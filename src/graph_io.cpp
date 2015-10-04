#include "graph_io.hpp"

using namespace seqan;


void printIdMap(boost::unordered_map<std::string, long> &my_id_map)
{
  // std::cout << my_id_map.size() << std::endl;
  std::string max_id;
  long max_count = 0;

  for (boost::unordered_map<std::string, long>::iterator it = my_id_map.begin() ; it != my_id_map.end() ; ++it)
  {
    printf("%8s: %4ld\n", it->first.c_str(), it->second);
    if (it->second > max_count)
    {
      max_count = it->second;
      max_id = it->first.c_str();
    }
  }
  std::cout << max_id << " wins!" << std::endl;
}


void writeDotFile(TGraph graph)
{
  std::ofstream dotFile("graph.dot");
  writeRecords(dotFile, graph, DotDrawing());

  dotFile.close();
}


void printGraph(TGraph graph)
{
  std::cout << graph << std::endl;
}


void
changeIdFormat(std::vector<std::string> &ids, std::string &gene, bool dont_change_id)
{
  // std::string asterisk();
  for (unsigned i = 0 ; i < ids.size() ; ++i)
  {
    std::size_t found = ids[i].find_first_of('*');
    gene = ids[i].substr(0,found);
    if (dont_change_id)
    {
      return;
    }
    ids[i] = ids[i].substr(found+1);

    // std::cout << ids[i].at(0) << std::endl;
    if (ids[i].at(0) == '0')
    {
      ids[i] = ids[i].substr(1);
    }

    while (true)
    {
      std::size_t found_colon = ids[i].find_first_of(':');
      // std::cout << found_colon << std::endl;
      if (found_colon == std::string::npos)
        break;

      ids[i] = ids[i].replace(ids[i].begin()+found_colon, ids[i].begin()+found_colon+1, "");
    }
    // std::cout << ids[i] << std::endl;
  }
}


void
getVcfHeader(VcfHeader &header)
{
  time_t rawtime;
  tm* timeinfo;
  char buffer [80];
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,10,"%Y%m%d",timeinfo);
  std::string date(buffer);
  appendValue(header, VcfHeaderRecord("fileformat", "VCFv4.1"));
  appendValue(header, VcfHeaderRecord("fileDate", date));
  appendValue(header, VcfHeaderRecord("source", "GyperV0.1"));
  appendValue(header, VcfHeaderRecord("reference", "genome.fa"));
  appendValue(header, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
  appendValue(header, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=String,Description=\"Phred-scaled genotype likelihoods\">"));
}


long adjustPositions(
  unsigned const & minSeqLen, // 10 - 90 (10)
  double const & beta,        // 0.1 - 0.9 (0.1)
  int const & bpQclip,        // 20 - 40 (5)
  int const & bpQskip         // 20 - 40 (5)
)
{
  long minSeq_adjust = (minSeqLen-10)/10; // 0, 1, 2, ..., 8
  long beta_adjust = (long) std::floor((beta*100.0-10.0)/10.0*9.0+0.5); // 0, 9, 18, ..., 72
  long bpQclip_adjust = (bpQclip - 20)/5*(72+8+1); // 81, ..., 324
  long bpQskip_adjust = (bpQskip - 20)/5*(324+72+8+1); // 405, ..., 
  std::cout << "minSeq_adjust = " << minSeq_adjust << std::endl;
  std::cout << "beta_adjust = " << beta_adjust << std::endl;
  std::cout << "bpQclip_adjust = " << bpQclip_adjust << std::endl;
  std::cout << "bpQskip_adjust = " << bpQskip_adjust << std::endl;
  std::cout << "adjust_positions = " << minSeq_adjust + beta_adjust + bpQclip_adjust + bpQskip_adjust << std::endl;
  return minSeq_adjust + beta_adjust + bpQclip_adjust + bpQskip_adjust;
}


long
getGeneBeginPos(
  std::string const & gene,
  unsigned const & minSeqLen,
  double const & beta,
  int const & bpQclip,
  int const & bpQskip
)
{
  int new_bpQclip = bpQclip;
  if (new_bpQclip == 37)
  {
    // Needed to avoid adjustments conflicts
    new_bpQclip = 40;
  }

  int new_bpQskip = bpQskip;
  if (new_bpQskip == 37)
  {
    // Needed to avoid adjustments conflicts
    new_bpQskip = 40;
  }

  if (gene == "HLA-A")
    return (long) 29942470 + adjustPositions(minSeqLen, beta, new_bpQclip, new_bpQskip); 
  if (gene == "HLA-B")
    return (long) 31353872 + adjustPositions(minSeqLen, beta, new_bpQclip, new_bpQskip);
  if (gene == "HLA-C")
    return (long) 31268749 + adjustPositions(minSeqLen, beta, new_bpQclip, new_bpQskip);
  if (gene == "HLA-DQA1")
    return (long) 32637406 + adjustPositions(minSeqLen, beta, new_bpQclip, new_bpQskip);
  if (gene == "HLA-DQB1")
    return (long) 32659464 + adjustPositions(minSeqLen, beta, new_bpQclip, new_bpQskip);
  if (gene == "HLA-DRB1")
    return (long) 32578770 + adjustPositions(minSeqLen, beta, new_bpQclip, new_bpQskip);

  std::cerr << "Warning: Gene " << gene << " not found!" << std::endl;
  return 0;
}

void
writeVcfFile(CharString & file_name,
             std::string const & pn,
             std::vector<std::string> ids,
             std::vector<std::vector<double> > & seq_scores,
             double const & max_score
            )
{
  writeVcfFile(file_name, pn, ids, seq_scores, max_score, 10, 0.1, 20, 20);
}

void
writeVcfFile(CharString & file_name,
             std::string const & pn,
             std::vector<std::string> ids,
             std::vector<std::vector<double> > & seq_scores,
             double const & max_score,
             unsigned const & minSeqLen,
             double const & beta,
             int const & bpQclip,
             int const & bpQskip
            )
{
  using namespace seqan;

  // std::ofstream stream_out(file_name, std::ofstream::out);
  std::ofstream myfile;
  myfile.open(toCString(file_name), std::ios_base::trunc); 
  // stream_out.open(file_name);

  if (ids.size() <= 2)
  {
    std::cerr << "Error: Size of ids is too small: " << ids.size() << "." << std::endl;
    return;
  }
  if (seq_scores.size() != ids.size())
  {
    std::cerr << "Error: The size of the seq_scores and ids should be the same." << std::endl;
    return;
  }
  if (pn.size() != 7)
  {
    std::cerr << "Error: PN should be 7 letters but not " << pn.size() << "." << std::endl;
    return;
  }

  std::ostringstream my_vcf_ss;
  VcfFileOut out(my_vcf_ss, Vcf());
  appendValue(contigNames(context(out)), "chr6");
  appendValue(sampleNames(context(out)), pn);

  VcfHeader header;
  getVcfHeader(header);
  writeHeader(out, header);

  std::string gene;
  // changeIdFormat(ids, gene, true);
  {
    std::size_t found = ids[0].find_first_of('*');
    gene = ids[0].substr(0,found);
  }

  VcfRecord record;
  record.rID = 0;
  record.beginPos = getGeneBeginPos(gene, minSeqLen, beta, bpQclip, bpQskip);
  record.id = gene;
  record.filter = "PASS";
  record.info = ".";
  record.format = "GT:PL";
  // record.ref = "0";
  record.ref = ids[0];

  // CharString alts(ids[1]);
  std::ostringstream my_alts;
  // my_alts << "1";
  my_alts << ids[1];
  

  for (unsigned i = 2 ; i < ids.size() ; ++i)
  {
    my_alts << ",";
    // my_alts << i;
    my_alts << ids[i];
  }

  CharString alts = my_alts.str();

  record.alt = alts;
  int i_max = -1;
  int j_max = -1;
  int phred_penalty = 20;
  std::ostringstream my_ss;

  for (unsigned i = 0 ; i < seq_scores.size() ; ++i)
  {
    for (unsigned j = 0 ; j <= i ; ++j)
    {
      if (i != 0 || j != 0)
      {
        my_ss << ",";
      }

      seq_scores[i][j] = max_score - seq_scores[i][j];
      if (i_max == -1 && seq_scores[i][j] < 1e-8)
      {
        i_max = i;
        j_max = j;
      }

      int phred_score = std::floor(seq_scores[i][j]*phred_penalty+0.5);
      if (phred_score > 255)
      {
        phred_score = 255;
      }
      char buffer [10];
      snprintf(buffer, 10, "%d", phred_score);
      CharString score(buffer);
      my_ss << score;
    } 
  }
  // CharString phred_scores = ":";
  // append(phred_scores, my_ss.str());

  char buffer [10];
  snprintf(buffer, 10, "%d/%d", j_max, i_max);
  CharString genoTypeInfoStr(buffer);
  // CharString genoTypeInfoStr = "";
  // append(genoTypeInfoStr, ids[j_max]);
  // append(genoTypeInfoStr, "/");
  // append(genoTypeInfoStr, ids[i_max]);
  append(genoTypeInfoStr, ":");
  append(genoTypeInfoStr, my_ss.str());
  appendValue(record.genotypeInfos, genoTypeInfoStr);
  writeRecord(out, record);

  myfile << my_vcf_ss.str();
  myfile.close();
}


void
testVcfFile(CharString file_name,
            std::string pn,
            std::vector<std::string> ids,
            std::vector<std::vector<double> > seq_scores,
            double max_score
            )
{
  writeVcfFile(file_name, pn, ids, seq_scores, max_score, 60, 0.6, 30, 25);
}

#include <stdio.h>
#include <ctime>
#include <fstream>
#include <string>

#include "graph_align.hpp"
#include "graph_builder.hpp"
#include "graph_io.hpp"
#include "graph.hpp"

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
// #include <seqan/bam_io/bam_index_bai.h>

#include <boost/algorithm/string/find.hpp>

struct callOptions
{
 public:
  CharString beta_list;
  std::vector<double> beta;
  int bpQclip;
  int bpQskip;
  CharString bamFile;
  CharString bam2;
  CharString bam3;
  CharString bam4;
  CharString gene;
  CharString minSeqLen_list;
  std::vector<unsigned> minSeqLen;
  // std::vector<int> minSeqs;
  CharString outputFolder;
  CharString vcfFile;
  CharString vcfOutputFolder;
  bool verbose;
  bool align_all_reads;
  bool thousand_genomes;
  int read_gap;

  void output(){
    std::cout << "CO.bamFile " << bamFile << std::endl;
    std::cout << "CO.bam2 " << bam2 << std::endl;
    std::cout << "CO.bam3 " << bam3 << std::endl;
    std::cout << "CO.bam4 " << bam4 << std::endl;
    std::cout << "CO.beta_list " << beta_list << std::endl;
    std::cout << "CO.beta.size() " << beta.size() << std::endl;
    std::cout << "CO.bpQclip " << bpQclip << std::endl;
    std::cout << "CO.bpQskip " << bpQskip << std::endl;
    std::cout << "CO.gene " << gene << std::endl;
    std::cout << "CO.minSeqLen_list " << minSeqLen_list << std::endl;
    std::cout << "CO.minSeqLen.size() " << minSeqLen.size() << std::endl;
    std::cout << "CO.outputFolder " << outputFolder << std::endl;
    std::cout << "CO.vcfFile " << vcfFile << std::endl;
    std::cout << "CO.vcfOutputFolder" << vcfOutputFolder << std::endl;
    std::cout << "CO.verbose " << verbose << std::endl;
    std::cout << "CO.align_all_reads " << align_all_reads << std::endl;
    std::cout << "CO.thousand_genomes " << thousand_genomes << std::endl;
    std::cout << "CO.read_gap " << read_gap << std::endl;
  };

callOptions():
  beta_list("0.6"), bpQclip(30), bpQskip(25), gene("DQA1"), minSeqLen_list("60"), outputFolder(), vcfOutputFolder(), verbose(false), align_all_reads(false), thousand_genomes(false), read_gap(1000) {}
};


template<typename T>
void parseArgList (CharString &args_in, std::vector<T> &vector_out)
{
  unsigned j = 0;
  for (unsigned i = 0 ; i < length(args_in) ; ++i)
  {
    if (args_in[i] != ',')
      continue;

    CharString copy(args_in);
    copy = infix(copy, j, i);
    std::string new_beta(toCString(copy));
    T new_beta_casted = boost::lexical_cast<T>(new_beta);
    j = i+1;
    vector_out.push_back(new_beta_casted);
  }

  CharString copy(args_in);
  copy = suffix(copy, j);
  std::string new_beta(toCString(copy));
  T new_beta_casted = boost::lexical_cast<T>(new_beta);
  vector_out.push_back(new_beta_casted);
}


ArgumentParser::ParseResult parseCommandLine(callOptions& CO, ArgumentParser & parser, int argc, char const ** argv )
{
  // setShortDescription(parser, "Converts a bam file into a fastq file.");
  setVersion(parser, "0.1");
  setDate(parser, "July 2015");

  // Main options
  addSection(parser, "Main options");
  addOption(parser, ArgParseOption("o", "output", "Outfile folder location ", ArgParseArgument::STRING, "output"));
  addOption(parser, ArgParseOption("vcf", "vcf_output", "Outfile folder location ", ArgParseArgument::STRING, "vcf_output"));
  addOption(parser, ArgParseOption("ms", "minSeqLen", "Minimum sequence length ", ArgParseArgument::STRING, "STRING"));
  addOption(parser, ArgParseOption("qc", "bpQclip", "Quality Clip Threshold for clipping. ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("qs", "bpQskip", "Quality Clip Threshold for skipping. ", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("v", "verbose", "Verbosity flag"));
  addOption(parser, ArgParseOption("a", "align_all_reads", "If used Gyper will use every read pair in the BAM file provided. Should only be used for very small BAMs, unless you have huge amount of RAM and patience."));
  addOption(parser, ArgParseOption("1k", "thousand_genomes", "Use reference from the 1000 Genomes project."));
  addOption(parser, ArgParseOption("b",
    "beta",
    "A scoring parameter for heterozygous scores. Higher beta means more likely to pick heterozugous.",
    ArgParseArgument::STRING, "STRING"));

  setDefaultValue(parser, "beta", "0.6");
  addOption(parser, ArgParseOption("b2", "bam2", "Second bam file", ArgParseArgument::INPUT_FILE, "bam2"));
  addOption(parser, ArgParseOption("b3", "bam3", "Third bam file", ArgParseArgument::INPUT_FILE, "bam3"));
  addOption(parser, ArgParseOption("b4", "bam4", "Fourth bam file", ArgParseArgument::INPUT_FILE, "bam4"));
  addOption(parser, ArgParseOption("rg", "read_gap", "The amount af bp to extend the read gap (i.e. the gap where reads are considered.) ", ArgParseArgument::INTEGER, "INT"));
  
  // setDefaultValue( parser, "qual", CO.bpQclip );

  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // Only extract  options if the program will continue after parseCommandLine()
  if (res != ArgumentParser::PARSE_OK)
    return res;

  if (isSet(parser, "minSeqLen"))
    getOptionValue(CO.minSeqLen_list, parser, "minSeqLen");
  if (isSet(parser, "output" ))
    getOptionValue(CO.outputFolder, parser, "output");
  if (isSet(parser, "vcf_output" ))
    getOptionValue(CO.vcfOutputFolder, parser, "vcf_output");
  if (isSet(parser, "bpQclip" ))
    getOptionValue(CO.bpQclip, parser, "bpQclip");
  if (isSet(parser, "bpQskip" ))
    getOptionValue(CO.bpQskip, parser, "bpQskip");
  if (isSet(parser, "bam2" ))
    getOptionValue(CO.bam2, parser, "bam2");
  if (isSet(parser, "bam3" ))
    getOptionValue(CO.bam3, parser, "bam3");
  if (isSet(parser, "bam4" ))
    getOptionValue(CO.bam4, parser, "bam4");
  if (isSet(parser, "beta"))
    getOptionValue(CO.beta_list, parser, "beta");
  if (isSet(parser, "align_all_reads"))
    getOptionValue(CO.align_all_reads, parser, "align_all_reads");
  if (isSet(parser, "thousand_genomes"))
    getOptionValue(CO.thousand_genomes, parser, "thousand_genomes");
  if (isSet(parser, "read_gap"))
    getOptionValue(CO.read_gap, parser, "read_gap");

  // Check if multiple values very given in a list
  parseArgList(CO.beta_list, CO.beta);
  parseArgList(CO.minSeqLen_list, CO.minSeqLen);
  getArgumentValue( CO.gene, parser, 0);
  getArgumentValue( CO.bamFile, parser, 1);
  CO.verbose = isSet(parser, "verbose");

  // Remove the HLA- part from genes if it is there
  if (CO.gene == "HLA-A")
  {
    CO.gene = "HLAA";
  }
  else if (CO.gene == "HLA-B")
  {
    CO.gene = "HLAB";
  }
  else if (CO.gene == "HLA-C")
  {
    CO.gene = "HLAC";
  }
  else if (CO.gene == "HLA-DQA1")
  {
    CO.gene = "DQA1";
  }
  else if (CO.gene == "HLA-DQB1")
  {
    CO.gene = "DQB1";
  }
  else if (CO.gene == "HLA-DRB1")
  {
    CO.gene = "DRB1";
  }

  if (CO.verbose)
  {
    CO.output();
  }

  return ArgumentParser::PARSE_OK;
}

typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

unsigned
getChromosomeRid(TBamContext const & bamContext, CharString const & chromosomeName)
{
  for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
  {
    if (contigNames(bamContext)[i] == chromosomeName)
    {
      return i;
    }
  }
  // std::cerr << "No chromosome '" << chromosomeName << "' found!" << std::endl;
  // std::cerr << "Available chromosome names are:\n";
  // for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
  // {
  //   std::cerr << contigNames(bamContext)[i] << std::endl;
  // }
  return 99999;
}

std::string
getRegion(callOptions & CO,
          std::string gene)
{
  std::ostringstream region_file;
  
  if (CO.thousand_genomes)
  {
    region_file << "data/region1k/" << gene;
  }
  else
  {
    region_file << "data/region/" << gene;
  }
  
  std::ifstream myfile(region_file.str());
  std::string line;

  if (myfile.is_open())
  {
    while (getline(myfile,line))
    {
      return line;
    }
    myfile.close();
  }
  std::cerr << "No region file found!" << std::endl;
  return line;
}


bool
validateRecord(callOptions & CO,
               TBamContext const & bamContext,
               BamAlignmentRecord const & record,
               std::string region
              )
{
  std::string delimiter = " ";
  size_t pos = 0;
  std::string token;
  while ((pos = region.find(delimiter)) != std::string::npos)
  {
    token = region.substr(0, pos);

    size_t colon_pos = token.find(":");
    std::string chromosomeName = token.substr(0, colon_pos);
    unsigned rID = getChromosomeRid(bamContext, chromosomeName);
    if (rID == 99999)
    {
      // This means we didn't find the chromosome
      rID = getChromosomeRid(bamContext, chromosomeName.substr(3));
      if (rID == 99999)
      {
        std::cerr << "I can't find the rID!" << std::endl;
      }
    }
    // std::cout << "chromosomeName = " << chromosomeName << " rID = " << rID << " record.rID = " << record.rID << std::endl;
    if (rID == (unsigned) record.rID)
    {
      size_t hyphen_pos = token.find("-");
      std::string begin_pos = token.substr(colon_pos + 1, hyphen_pos - (colon_pos + 1));
      std::string end_pos = token.substr(hyphen_pos + 1);
      if ((unsigned) record.beginPos > stoul(begin_pos) - length(record.seq) - CO.read_gap && (unsigned) record.beginPos < stoul(end_pos) + CO.read_gap)
      {
        return true;
      }
      // else
      // {
      //   std::cout << "NOPE: record.beginPos = " << record.beginPos << " begin_pos = " << begin_pos << " end_pos = " << end_pos << std::endl;
      // }
    }

    // std::cout << token << std::endl;
    region.erase(0, pos + delimiter.length());
  }
  // std::cout << "NOPE: record.rID = " << record.rID << " record.beginPos = " << record.beginPos << std::endl;
  return false;
}


bool
jumpAround(callOptions & CO,
           BamFileIn & bamFileIn,
           TBamContext const & bamContext,
           BamIndex<Bai> const & index,
           std::string & region,
           unsigned & begin_pos,
           unsigned & end_pos
          )
{
  if (region.size() < 6)
    return false;

  std::string delimiter = " ";
  size_t pos = region.find(delimiter);
  std::string token;
  token = region.substr(0, pos);
  // std::cout << "token = " << token << std::endl;

  size_t colon_pos = token.find(":");
  std::string chromosomeName = token.substr(0, colon_pos);
  unsigned rID = getChromosomeRid(bamContext, chromosomeName);
  
  if (rID == 99999)
  {
    // This means we didn't find the chromosome as chrX, so let's try if we can find just X
    rID = getChromosomeRid(bamContext, chromosomeName.substr(3));
    if (rID == 99999)
    {
      std::cerr << "I can't find the rID!" << std::endl;
    }
  }

  size_t hyphen_pos = token.find("-");
  std::string begin_str = token.substr(colon_pos + 1, hyphen_pos - (colon_pos + 1));
  // std::cout << "begin_str = " << begin_str << std::endl;

  begin_pos = stoul(begin_str);
  std::string end_str = token.substr(hyphen_pos + 1);
  // std::cout << "end_str = " << end_str << std::endl;
  end_pos = stoul(end_str);
  region.erase(0, pos + delimiter.length());

  bool hasAlignments;

  if (!jumpToRegion(bamFileIn, hasAlignments, rID, begin_pos, end_pos, index))
  {
    std::cerr << "Couldn't jump in bam file!" << std::endl;
    return false;
  }
  else if (CO.verbose)
  {
    std::cout << "Jumping to positions " << chromosomeName << ":" << begin_pos << "-" << end_pos << " in the bam file." << std::endl;
  }


  if (!hasAlignments)
  {
    std::cerr << "No alignments found at that location!" << std::endl;
    begin_pos = 0;
    end_pos = 0;
  }

  return true;
}


int
addToBars(callOptions & CO,
          std::map<CharString, BamAlignmentRecord> & bars1,
          std::map<CharString, BamAlignmentRecord> & bars2,
          CharString & bam,
          std::string & region
         )
{
  if (bam == "")
    return 0;

  if (CO.verbose)
  {
    std::cout << "Adding " << bam << std::endl;
  }
  
  std::string region_copy(region);
  char* bam_file_name = toCString(bam);
  BamFileIn bamFileIn(bam_file_name);
  BamHeader header;
  readHeader(header, bamFileIn);
  
  TBamContext const & bamContext = context(bamFileIn);
  // std::cout << "Chromosome 6 has rID: " << getChromosomeRid(bamContext, "chr6") << std::endl;
  
  BamAlignmentRecord record;

  if (!CO.align_all_reads)
  {
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, strcat(bam_file_name, ".bai")))
    {
      std::cerr << "ERROR: Could not read BAI index file " << bam_file_name << "\n";
      return 1;
    }

    while (true)
    {
      unsigned begin_pos;
      unsigned end_pos;

      if(!jumpAround(CO, bamFileIn, bamContext, baiIndex, region, begin_pos, end_pos))
      {
        break;
      }
      
      readRecord(record, bamFileIn);

      while((unsigned) record.beginPos < begin_pos - CO.read_gap - length(record.seq))
      {
        readRecord(record, bamFileIn);
      }

      while((unsigned) record.beginPos < end_pos + CO.read_gap && !atEnd(bamFileIn))
      {
        // std::cout << "record.beginPos = " << record.beginPos << std::endl;
        if(!validateRecord(CO, bamContext, record, region_copy))
        {
          readRecord(record, bamFileIn);
          continue;
        }

        if (hasFlagFirst(record))
        {
          if (bars1.count(record.qName) == 1)
          {
            readRecord(record, bamFileIn);
            continue;
          }

          bars1[record.qName] = record;
        }
        else
        {
          if (bars2.count(record.qName) == 1)
          {
            readRecord(record, bamFileIn);
            continue;
          }

          bars2[record.qName] = record;
        }

        readRecord(record, bamFileIn);
      }
    }
  }
  else
  {
    while (!atEnd(bamFileIn))
    {
      readRecord(record, bamFileIn);

      // std::cout << "record.beginPos = " << record.beginPos;
      // std::cout << ": record.rID = " << record.rID << std::endl;

      if (hasFlagFirst(record))
      {
        if (bars1.count(record.qName) == 1)
        {
          continue;
        }

        bars1[record.qName] = record;
      }
      else
      {
        if (bars2.count(record.qName) == 1)
        {
          continue;
        }

        bars2[record.qName] = record;
      }
    }
  }
  
  return 0;
}


inline int
qualToInt( char c )
{
  return (int) c - 33;
}


boost::dynamic_bitset<>
trimReadEnds (DnaString & seq,
              CharString & qual,
              int const & bpQclip,
              int const & bpQskip
             )
{
  using namespace seqan;

  int beg, end;

  for (beg = 0; (beg < (int) length(seq)) && (qualToInt(qual[beg]) < bpQclip); ++beg);
  for (end = (int) length(seq) - 1; (end > 0) && (qualToInt(qual[end]) < bpQclip); --end);
  
  if (beg > end)
  {
    seq = "";
    qual = "";
  }
  else
  {
    seq = infix(seq, beg, end+1);
    qual = infix(qual, beg, end+1);
  }

  boost::dynamic_bitset<> quality_bitset(length(seq));

  for (unsigned pos = 0 ; pos < length(seq) ; ++pos)
  {
    if (qualToInt(qual[pos]) < bpQskip)
    {
      quality_bitset[pos] = 1;
    }
  }
  return quality_bitset;
}


void
appendToFile (std::stringstream& ss, CharString &fileName)
{  
  using namespace std;
  string myString = ss.str();  
  ofstream myfile;  
  myfile.open(toCString(fileName), ios_base::app); 
  myfile << myString;  
  myfile.close();  
}


void
writeToFile (std::stringstream& ss, CharString &fileName)
{
  using namespace std;
  string myString = ss.str();
  ofstream myfile;
  myfile.open(toCString(fileName), ios_base::trunc); 
  myfile << myString;
  myfile.close();
}


std::string
fourDigit(callOptions & CO, std::string my_id)
{
  size_t n = std::count(my_id.begin(), my_id.end(), ':');

  unsigned colons;
  // if (CO.gene == "HLAB" || CO.gene == "HLAC")
  // {
  //   // Only 2 digit resolution for HLAB and HLAC
  //   colons = 0;
  // }
  // else
  // {
  //   colons = 1;
  // }
  colons = 1;

  if (n > colons)
  {
    boost::iterator_range<std::string::iterator> r = boost::find_nth(my_id, ":", colons);
    size_t distance = std::distance(my_id.begin(), r.begin());
    my_id = my_id.substr(0, distance);
  }

  char last_char = my_id.back();

  if (!isdigit(last_char))
  {
    std::string no_char_id(my_id.begin(), my_id.end() - 1);
    return no_char_id;
  }

  return my_id;
}


size_t
findAvailableGenotypes(callOptions & CO, boost::unordered_set<std::string> & available_alleles)
{
  std::ostringstream file_location;

  if (CO.thousand_genomes)
  {
    file_location << "data/available_alleles/1000genomes/";
  }
  else
  {
    file_location << "data/available_alleles/icelandic_alleles/";
  }

  file_location << CO.gene;
  std::ifstream available_genotypes_file_stream(file_location.str());
  std::string line;

  available_genotypes_file_stream >> line;
  size_t n = std::count(line.begin(), line.end(), ':');
  available_alleles.insert(line);

  while (available_genotypes_file_stream >> line)
  {
    available_alleles.insert(line);
  }

  available_genotypes_file_stream.close();
  return n;
}


void
handleOutput (callOptions & CO,
              const double & beta,
              // const unsigned & min_seq_index,
              const std::string & pn,
              std::vector<std::string> & ids,
              std::vector<std::vector<double> > seq_scores,
              double const & total_matches
              )
{
  // Four digit vector
  std::vector<std::string> ids_four;

  // The 4 digit index numbers
  boost::unordered_map<std::string, unsigned> four;
  unsigned j = 0;

  for (unsigned i = 0; i < ids.size(); ++i)
  {
    // std::cout << "Id before: " << my_id << std::endl;
    std::string my_id = fourDigit(CO, ids[i]);
    // std::cout << "Id after: " << my_id << std::endl;

    if (!four.count(my_id))
    {
      four[my_id] = j;
      ids_four.push_back(my_id);
      ++j;
    }
  }

  std::vector<std::vector<double> > seq_scores_four;
  boost::unordered_map<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > index_map;

  {
    std::vector<double> inner_seq_scores;
    inner_seq_scores.resize(four.size());
    std::vector<std::vector<double> > outer_seq_scores;
    seq_scores_four.resize(four.size(), inner_seq_scores);
  }

  boost::unordered_set<std::string> available_alleles;
  unsigned available_digits = findAvailableGenotypes(CO, available_alleles);
  std::cout << "available_digits = " << available_digits << std::endl;

  std::string short_id_i;
  std::string short_id_j;

  for (unsigned i = 0 ; i < ids.size() ; ++i)
  {
    for (unsigned j = 0 ; j <= i ; ++j)
    {
      // No need to scale the homozygous solutions
      if (i != j)
      {
        seq_scores[i][j] *= beta;
        seq_scores[i][j] += (seq_scores[i][i]+seq_scores[j][j])/2.0*(1.0-beta);
      }

      {
        // Get the short ID for i and skip this iteration if it isn't available
        boost::iterator_range<std::string::iterator> r = boost::find_nth(ids[i], ":", available_digits);
        size_t distance = std::distance(ids[i].begin(), r.begin());
        short_id_i = ids[i].substr(0, distance);
        // std::cout << "short_id_i = " << short_id_i << std::endl;
      }

      if (available_alleles.find(short_id_i) == available_alleles.end())
      {
        continue;
      }

      {
        // Same for j
        boost::iterator_range<std::string::iterator> r = boost::find_nth(ids[j], ":", available_digits);
        size_t distance = std::distance(ids[j].begin(), r.begin());
        short_id_j = ids[j].substr(0, distance);
        // std::cout << "short_id_j = " << short_id_j << std::endl;
      }

      if (available_alleles.find(short_id_j) == available_alleles.end())
      {
        continue;
      }

      unsigned i_four = four[fourDigit(CO, ids[i])];
      unsigned j_four = four[fourDigit(CO, ids[j])];

      if (seq_scores[i][j] > seq_scores_four[i_four][j_four])
      {
        seq_scores_four[i_four][j_four] = seq_scores[i][j];
        std::pair<unsigned, unsigned> index_pair(i, j);
        std::pair<unsigned, unsigned> index_pair_four(i_four, j_four);
        // std::cout << "Creating a index pair with: " << ids[i] << ", " << ids[j] << std::endl;
        index_map[index_pair_four] = index_pair;
      }
    }
  }

  unsigned i_max = 0;
  unsigned j_max = 0;
  double max_score = 0;
  
  for (unsigned i = 0 ; i < ids.size() ; ++i)
  {
    // std::cout << std::setw(17) << ids[i] << ": ";
    for (unsigned j = 0 ; j < ids.size() ; ++j)
    {
      // printf("%7.2f ", seq_scores[i][j]);
      
      if (seq_scores[i][j] > max_score)
      {
        i_max = i;
        j_max = j;
        max_score = seq_scores[i][j];
      }
    }
    // std::cout << std::endl;
  }
  
  unsigned i_max_short = 0;
  unsigned j_max_short = 0;
  double max_score_short = 0;

  // Only print this if it's not absurdly large (and in verbose mode)
  for (unsigned i = 0; i < ids_four.size(); ++i)
  {
    std::cout << std::setw(15) << ids_four[i] << ": ";
    for (unsigned j = 0; j < ids_four.size(); ++j)
    {
      printf("%6.1f ", seq_scores_four[i][j]);

      if (seq_scores_four[i][j] > max_score_short)
      {
        i_max_short = i;
        j_max_short = j;
        max_score_short = seq_scores_four[i][j];
      }
    }
    std::cout << std::endl;
  }
  
  // Get the original alleles
  std::pair<unsigned, unsigned> best_alleles(i_max_short, j_max_short);
  std::pair<unsigned, unsigned> allele_ids = index_map[best_alleles];

  if (CO.verbose)
  {
    std::cout << "Total matches are: " << total_matches << std::endl;
  }

  if (total_matches != 0)
  {
    std::cout << "Oracle: Out of all the alleles, I'd say this person has   " << ids[i_max] << " and " << ids[j_max] << " with a score " << max_score << "!" << std::endl;
    std::cout << "Oracle: Out of available alleles, I'd say this person has " << ids[allele_ids.first] << " and " << ids[allele_ids.second] << " with a score " << max_score_short << "!" << std::endl;
  
    // if (i_max != j_max)
    // {
    //   std::cout << " (pre-adjustment: " << (max_score-(seq_scores[i_max][i_max]+seq_scores[j_max][j_max])*(1-beta)/2.0)/beta << " or ";
    //   std::cout << 100.0*(max_score-(seq_scores[i_max][i_max]+seq_scores[j_max][j_max])*(1-beta)/2.0)/beta/total_matches << "%)";
    // }
    // else
    // {
    //   std::cout << " (" << 100.0*max_score/total_matches << "%)";
    // }
    // std::cout << std::endl;
    
    {
      CharString outputFolder;
      if (CO.outputFolder == "")
      {
        outputFolder = "output/gyper/";
      }
      else
      {
        outputFolder = CO.outputFolder;
      }

      append(outputFolder, CO.gene);
      append(outputFolder, "/");
      append(outputFolder, pn);
      append(outputFolder, ".txt");

      std::stringstream my_ss;
      my_ss << ids[allele_ids.first] << "\t" << ids[allele_ids.second];
      // my_ss << ids[i_max] << "\t" << ids[j_max];
      // my_ss << ids[i_max] << "\t" << ids[j_max] << "\t" << CO.minSeqLen[min_seq_index];
      // my_ss << "\t" << beta << "\t" << CO.bpQclip << "\t" << CO.bpQskip << std::endl;
      std::cout << "Gyper output: " << outputFolder << std::endl;
      writeToFile(my_ss, outputFolder);
    }
  }

  // Print VCF
  {
    if (CO.verbose)
    {
      // if (CO.gene == "HLAB" || CO.gene == "HLAC")
      // {
      //   std::cout << "Number of 2 digit ids are: " << ids_four.size() << std::endl;
      // }
      // else
      // {
        std::cout << "Number of 4 digit ids are: " << ids_four.size() << std::endl;
      // }
    }

    CharString vcfOutputFolder;
    if (CO.vcfOutputFolder == "")
    {
      vcfOutputFolder = "output/VCF/";
    }
    else
    {
      vcfOutputFolder = CO.vcfOutputFolder;
    }
  
    append(vcfOutputFolder, CO.gene);
    append(vcfOutputFolder, "/");
    // append(vcfOutputFolder, boost::lexical_cast<std::string>(CO.minSeqLen[min_seq_index]));
    // append(vcfOutputFolder, "/");
    // std::ostringstream ss;
    // ss << std::fixed << std::setprecision(2) << beta;
    // append(vcfOutputFolder, ss.str());
    // append(vcfOutputFolder, "/");
    // append(vcfOutputFolder, boost::lexical_cast<std::string>(CO.bpQclip));
    // append(vcfOutputFolder, "/");
    // append(vcfOutputFolder, boost::lexical_cast<std::string>(CO.bpQskip));
    // append(vcfOutputFolder, "/");
    append(vcfOutputFolder, pn);
    append(vcfOutputFolder, ".vcf");
    std::cout << "VCF output: " << vcfOutputFolder << std::endl;
    // writeVcfFile(vcfOutputFolder,
    //              pn, ids,
    //              seq_scores,
    //              max_score,
    //              CO.minSeqLen[min_seq_index],
    //              beta,
    //              CO.bpQclip,
    //              CO.bpQskip
    //             );
    writeVcfFile(vcfOutputFolder,
                 pn,
                 ids_four,
                 seq_scores_four,
                 max_score_short
                );
  }
}


int main (int argc, char const ** argv)   
{
  // Let's go!
  clock_t begin = clock();

  callOptions CO;

  ArgumentParser parser("gyper");
  addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIGENE\\fP\" \"\\fIBAMFILE\\fP\"");
  addDescription(parser, "Graph genotYPER. Maps sequencing data to possible alleles using a partial order graph.");

  // Requires at least one bam file as argument.
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "gene"));
  addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "bamfile"));

  if (parseCommandLine(CO, parser, argc, argv) != ArgumentParser::PARSE_OK)
    return ArgumentParser::PARSE_ERROR;

  std::string bam_file_name_str(toCString(CO.bamFile));
  std::string forward_slash("/");
  std::string pn;

  if (bam_file_name_str.find(forward_slash) != std::string::npos)
  {
    std::string::iterator str_it = std::find_end(bam_file_name_str.begin(), bam_file_name_str.end(), forward_slash.begin(), forward_slash.end())+1;
    char buffer[7];
    std::size_t length = bam_file_name_str.copy(buffer,7,str_it-bam_file_name_str.begin());
    buffer[length]='\0';
    pn.assign(buffer);
  }
  else
  {
    pn = bam_file_name_str.substr(0, 7);
  }

  std::string region = getRegion(CO, toCString(CO.gene));

  if (CO.verbose)
  {
    std::cout << "region = " << region << std::endl;
  }

  std::map< CharString, BamAlignmentRecord > bars1;
  std::map< CharString, BamAlignmentRecord > bars2;
 
  if (addToBars(CO, bars1, bars2, CO.bamFile, region))
    return 1;

  if (addToBars(CO, bars1, bars2, CO.bam2, region))
    return 1;

  if (addToBars(CO, bars1, bars2, CO.bam3, region))
    return 1;

  if (addToBars(CO, bars1, bars2, CO.bam4, region))
    return 1;

  if (CO.verbose)
  {
    std::cout << "Bam alignment records:" << std::endl;
    std::cout << "Before: bars1.size() = " << bars1.size();
    std::cout << ", bars2.size() = " << bars2.size() << std::endl;
  }
  
  std::vector<CharString> no_mates;

  for (std::map<CharString, BamAlignmentRecord>::iterator it = bars1.begin() ; it != bars1.end() ; ++it)
  {
    if (bars2.count(it->first) == 0)
    {
      no_mates.push_back(it->first);
    }
  }

  for (std::map<CharString, BamAlignmentRecord>::iterator it = bars2.begin() ; it != bars2.end() ; ++it)
  {
    if (bars1.count(it->first) == 0)
    {
      no_mates.push_back(it->first);
    }
  }

  for (std::vector<CharString>::iterator it = no_mates.begin() ; it != no_mates.end() ; ++it)
  {
    bars1.erase(*it);
    bars2.erase(*it);
  }

  if (CO.verbose)
  {
    std::cout << "After : bars1.size() = " << bars1.size();
    std::cout << ", bars2.size() = " << bars2.size() << std::endl;
  }

  printf("[%6.2f] Read pairs filtered.\n", double(clock()-begin) / CLOCKS_PER_SEC);

  std::vector<std::string> ids;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;
  String<TVertexDescriptor> order;
  std::vector<ExactBacktracker> backtracker1;
  std::vector<ExactBacktracker> reverse_backtracker1;
  std::vector<ExactBacktracker> backtracker2;
  std::vector<ExactBacktracker> reverse_backtracker2;
  boost::unordered_set<TVertexDescriptor> free_nodes;
  free_nodes.insert(1);
  TGraph graph;
  std::string gene(toCString(CO.gene));

  if (gene == "DQA1")
  {
    createDqa1Graph(graph, vertex_vector, ids, edge_ids, free_nodes, order);
  }
  else if (gene == "DQB1")
  {
    createDqb1Graph(graph, vertex_vector, ids, edge_ids, free_nodes, order);
  }
  else if (gene == "DRB1")
  {
    createDrb1Graph(graph, vertex_vector, ids, edge_ids, free_nodes, order);
  }
  else if (gene == "HLAA")
  {
    createHlaaGraph(graph, vertex_vector, ids, edge_ids, free_nodes, order);
  }
  else if (gene == "HLAB")
  {
    createHlabGraph(graph, vertex_vector, ids, edge_ids, free_nodes, order);
  }
  else if (gene == "HLAC")
  {
    createHlacGraph(graph, vertex_vector, ids, edge_ids, free_nodes, order);
  }
  else
  {
    std::cerr << "I don't support gene " << gene << std::endl;
    return 1;
  }

  printf("[%6.2f] Graph created.\n", double(clock()-begin) / CLOCKS_PER_SEC);
  DnaString sequence1;
  DnaString sequence2;
  
  std::vector<std::vector<std::vector<double> > > seq_scores;
  
  {
    std::vector<double> inner_seq_scores;
    inner_seq_scores.resize(ids.size());
    std::vector<std::vector<double> > outer_seq_scores;
    outer_seq_scores.resize(ids.size(), inner_seq_scores);
    seq_scores.resize(CO.minSeqLen.size(), outer_seq_scores);
  }


  std::vector<double> total_matches;
  total_matches.resize(CO.minSeqLen.size(), 0.0);
  boost::unordered_map<std::string, int> highest_interval_map;
  boost::unordered_map<std::string, int> medium_interval_map;

  for (std::map<CharString, BamAlignmentRecord>::iterator it = bars1.begin() ; it != bars1.end() ; ++it)
  {
    sequence1 = it->second.seq;
    sequence2 = bars2[it->first].seq;

    boost::dynamic_bitset<> qual1 = trimReadEnds(sequence1, it->second.qual, CO.bpQclip, CO.bpQskip);
    boost::dynamic_bitset<> qual2 = trimReadEnds(sequence2, bars2[it->first].qual, CO.bpQclip, CO.bpQskip);
    std::vector<unsigned> seq_scores_indices;

    for (unsigned min_seq = 0; min_seq < CO.minSeqLen.size(); ++min_seq)
    {
      if (seqan::length(sequence1) >= CO.minSeqLen[min_seq] && seqan::length(sequence2) >= CO.minSeqLen[min_seq])
      {
        seq_scores_indices.push_back(min_seq);
      }
    }

    if (seq_scores_indices.size() == 0)
      continue;

    // if (CO.verbose)
    // {
    //   std::cout << "sequence1 = " << sequence1 << std::endl;
    //   std::cout << "sequence2 = " << sequence2 << std::endl;
    //   std::cout << "qual1 = " << it->second.qual << std::endl;
    //   std::cout << "qual2 = " << bars2[it->first].qual << std::endl;
    // }

    std::vector<TVertexDescriptor> matching_vertices1;
    std::vector<TVertexDescriptor> reverse_matching_vertices1;
    std::vector<TVertexDescriptor> matching_vertices2;
    std::vector<TVertexDescriptor> reverse_matching_vertices2;

    align_sequence(sequence1,
                   qual1,
                   graph,
                   vertex_vector,
                   order,
                   backtracker1,
                   reverse_backtracker1,
                   free_nodes,
                   matching_vertices1,
                   reverse_matching_vertices1
                  );

    if (matching_vertices1.size() == 0 && reverse_matching_vertices1.size() == 0)
    {
      continue;
    }

    align_sequence(sequence2,
                   qual2,
                   graph,
                   vertex_vector,
                   order,
                   backtracker2,
                   reverse_backtracker2,
                   free_nodes,
                   matching_vertices2,
                   reverse_matching_vertices2
                  );

    if (matching_vertices2.size() == 0 && reverse_matching_vertices2.size() == 0)
    {
      continue;
    }

    std::string read_group = toCString(myExtractTagValue(it->second.tags));

    // TODO: Add this as options for the program
    int highest_distance = 800;
    int medium_distance = 350;

    int best_score = highest_distance;
    TVertexDescriptor best_vertex1 = -1;
    TVertexDescriptor best_vertex2 = -1;

    bool reverse = false;
    for (std::vector<TVertexDescriptor>::iterator match_it = matching_vertices1.begin() ; match_it != matching_vertices1.end() ; ++match_it)
    {
      for (std::vector<TVertexDescriptor>::iterator match_it2 = matching_vertices2.begin() ; match_it2 != matching_vertices2.end() ; ++match_it2)
      {
        // std::cout << vertex_vector[*match_it].level << " " << vertex_vector[*match_it2].level << std::endl;
        int difference = abs(vertex_vector[*match_it].level - vertex_vector[*match_it2].level);
        if (difference < highest_distance)
        {
          int offset = abs(medium_distance - difference);
          best_score = offset;
          best_vertex1 = *match_it;
          best_vertex2 = *match_it2;
        }
      }
    }

    for (std::vector<TVertexDescriptor>::iterator match_it = reverse_matching_vertices1.begin() ; match_it != reverse_matching_vertices1.end() ; ++match_it)
    {
      for (std::vector<TVertexDescriptor>::iterator match_it2 = reverse_matching_vertices2.begin() ; match_it2 != reverse_matching_vertices2.end() ; ++match_it2)
      {
        // std::cout << vertex_vector[*match_it].level << " " << vertex_vector[*match_it2].level << std::endl;
        int difference = abs(vertex_vector[*match_it].level - vertex_vector[*match_it2].level);
        if (difference < highest_distance)
        {
          int offset = abs(medium_distance - difference);
          if (offset < best_score)
          {
            best_score = offset;
            best_vertex1 = *match_it;
            best_vertex2 = *match_it2;
            reverse = true;
          }
        }
      }
    }
    
    if (reverse)
    {
      backtracker1 = reverse_backtracker1;
      backtracker2 = reverse_backtracker2;
    }

    if (best_vertex1 == (TVertexDescriptor) -1)
    {
      continue;
    }

    boost::dynamic_bitset<> ids_found1 = backTrackAndCount(ids, backtracker1, best_vertex1, edge_ids);
    // std::cout << "ids_found1 = " << ids_found1 << std::endl;
    if (ids_found1.find_first() == ids_found1.npos)
    {
      // No bits are set in ids_found1
      continue;
    }

    boost::dynamic_bitset<> ids_found2 = backTrackAndCount(ids, backtracker2, best_vertex2, edge_ids);
    // std::cout << "ids_found2 = " << ids_found2 << std::endl;
    if (ids_found2.find_first() == ids_found1.npos)
    {
      // No bits are set in ids_found2
      continue;
    }

    boost::dynamic_bitset<> ids_intersection(ids_found1.size());
    ids_intersection = ids_found1 & ids_found2;
    std::vector<bool> hit;
    hit.resize(CO.minSeqLen.size(), false);

    for (std::vector<unsigned>::iterator seq_it = seq_scores_indices.begin(); seq_it != seq_scores_indices.end(); ++seq_it)
    {
      // Full reference alleles
      for (unsigned i = 0 ; i < ids.size() ; ++i)
      {
        for (unsigned j = 0 ; j <= i ; ++j)
        {
          if (ids_intersection[i] || ids_intersection[j])
          {
            seq_scores[*seq_it][i][j] += 1.0;
            hit[*seq_it] = true;
          }
        }
      }
      if (hit[*seq_it])
      {
        total_matches[*seq_it] += 1.0;
      }
    }    
  }

  printf("[%6.2f] Sequences aligned.\n", double(clock()-begin) / CLOCKS_PER_SEC);
  for (unsigned min_seq_index = 0; min_seq_index < CO.minSeqLen.size(); ++min_seq_index)
  {
    for (std::vector<double>::iterator it = CO.beta.begin() ; it != CO.beta.end() ; ++it)
    {
      printf("[%6.2f] Output with beta = %.2f and minSeqLen = %u\n", 
        double(clock()-begin) / CLOCKS_PER_SEC,
        *it,
        CO.minSeqLen[min_seq_index]
      );
      handleOutput(CO, *it, pn, ids, seq_scores[min_seq_index], total_matches[min_seq_index]);
    }
  }
  
  printf("[%6.2f] Done.\n", double(clock()-begin) / CLOCKS_PER_SEC);
  return 0;
}

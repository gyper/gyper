#include <gyper/constructor.hpp>

namespace gyper
{

Constructor::Constructor ()
{
	TGraph graph();
}

Constructor::Constructor (Options & CO)
{
  TGraph graph();
  Constructor::Constructor::CO = CO; // Replaces default options
}

void
Constructor::add_reference_sequence_preceding_a_point(TVertexDescriptor prev_vertex, unsigned point)
{
  bool sequence_started = false;
  unsigned starting_pos = vertex_labels[prev_vertex].order;
  seqan::String<seqan::Dna> current_dna = "";

  if (point > length(reference_sequence))
  {
    point = length(reference_sequence);
  }

  for (unsigned pos = starting_pos; pos < point; ++pos)
  {
    if (reference_sequence[pos] == 'N')
    {
      if (sequence_started)
      {
        TVertexDescriptor target_vertex = seqan::addVertex(graph);
        VertexLabel new_vertex_label(starting_pos, current_dna);
        seqan::addEdge(graph, prev_vertex, target_vertex);
        vertex_label_map[new_vertex_label] = target_vertex;
        vertex_labels.push_back(new_vertex_label);
        // std::cout << length(current_dna) << std::endl;
        current_dna = "";
        prev_vertex = target_vertex;
        sequence_started = false;
      }

      continue;
    }

    if (!sequence_started)
    {
      TVertexDescriptor target_vertex = seqan::addVertex(graph);
      VertexLabel new_vertex_label(pos, "");
      seqan::addEdge(graph, prev_vertex, target_vertex);
      vertex_label_map[new_vertex_label] = target_vertex;
      vertex_labels.push_back(new_vertex_label);
      prev_vertex = target_vertex;
      starting_pos = pos;
      sequence_started = true;
    }
    
    appendValue(current_dna, reference_sequence[pos]);
  }
}

void
Constructor::add_reference_sequence_to_graph(seqan::String<seqan::Dna5> & sequence)
{
  TVertexDescriptor prev_vertex = seqan::addVertex(graph);

  seqan::String<seqan::Dna> initial_dna = "";
  VertexLabel initial_vertex(0, initial_dna);

  vertex_label_map[initial_vertex] = prev_vertex;
  vertex_labels.push_back(initial_vertex);
  
  bool sequence_started = false;
  seqan::String<seqan::Dna> current_dna = "";
  unsigned starting_pos = 0;

  for (unsigned pos = 0; pos < length(sequence); ++pos)
  {
    if (sequence[pos] == 'N')
    {
      if (sequence_started)
      {
        TVertexDescriptor target_vertex = seqan::addVertex(graph);
        VertexLabel new_vertex_label(starting_pos, current_dna);
        seqan::addEdge(graph, prev_vertex, target_vertex);
        vertex_label_map[new_vertex_label] = target_vertex;
        vertex_labels.push_back(new_vertex_label);
        // std::cout << length(current_dna) << std::endl;
        current_dna = "";
        prev_vertex = target_vertex;
        sequence_started = false;
      }

      continue;
    }

    if (!sequence_started)
    {
      TVertexDescriptor target_vertex = seqan::addVertex(graph);
      VertexLabel new_vertex_label(pos, "");
      seqan::addEdge(graph, prev_vertex, target_vertex);
      vertex_label_map[new_vertex_label] = target_vertex;
      vertex_labels.push_back(new_vertex_label);
      prev_vertex = target_vertex;
      starting_pos = pos;
      sequence_started = true;
    }
    
    appendValue(current_dna, sequence[pos]);
  }
}

// void
// Gyper::create_reference_graph(seqan::String<char> region)
// {
//   // Extract sequence from a region in a FASTA file
//   seqan::String<seqan::Dna5> sequence = reference_fasta_ptr->extract_region(region);
//   add_reference_sequence_to_graph(sequence);
//   return;
// }

void
Constructor::create_HLA_graph()
{
  unsigned number_of_exons = get_number_of_exons();

  // Add p3 region
  add_FASTA_region(false, 0, false, true, false);

  // First exon
  {
    add_FASTA_region(true, number_of_exons, false, false, false);
    --number_of_exons;
  }

  while (number_of_exons >= 1)
  {
    // Intron
    add_FASTA_region(false, number_of_exons, true, false, false);

    // Exon
    if (CO.exon_2_and_3 && (number_of_exons != 2 && number_of_exons != 3 && number_of_exons != 4))
    {
      add_FASTA_region(false, number_of_exons, false, false, false);
    }
    else
    {
      add_FASTA_region(true, number_of_exons, false, false, false);
    }

    --number_of_exons;
  }

  // Final UTR
  add_FASTA_region(false, 0, false, false, true);

  // Finally, sort the nodes in topological order
  seqan::topologicalSort(order, graph);
}

void
Constructor::index()
{
  TKmerMap kmer_map;

  // for (Iterator<String<TVertexDescriptor const> const>::Type it = begin(order) ; it != end(order) ; ++it)
  // {
  //   TVertexDescriptor const & source_vertex = *it;

  //   if (free_nodes.count(source_vertex) == 0)
  //   {
  //     boost::dynamic_bitset<> id_bits(edge_ids.begin()->second.size());
  //     id_bits.flip();
  //     checkKmers(vertex_vector[source_vertex].dna, source_vertex, source_vertex, graph, vertex_vector, free_nodes, edge_ids, id_bits, kmer_map, static_cast<std::size_t>(CO.k));
  //   }
  // }

  // return kmer_map;
}

unsigned
Constructor::get_number_of_exons()
{
  // Set the number of exons
  if (CO.gene == "HLAA" || CO.gene == "HLAC")
  {
    return 8;
  }
  else if (CO.gene == "HLAB")
  {
    return 7;
  }
  else if (CO.gene == "DQA1")
  {
    return 4;
  }
  else if (CO.gene == "DQB1" || CO.gene == "DRB1")
  {
    return 6;
  }
  std::cerr << "Unsupported HLA gene: " << CO.gene << std::endl;
  return 0;
}

std::string
Constructor::get_HLA_base_path()
{
  std::stringstream base_path;
  base_path << gyper_SOURCE_DIRECTORY << "/data/haplotypes/hla/references/" << CO.gene << "/";
  return base_path.str();
}

void
Constructor::add_FASTA_region(bool add_bitstrings, int feature_number, bool intron_region, bool p3_region, bool p5_region)
{
  std::string base_path = get_HLA_base_path();

  if (p3_region)
  {
    base_path.append("p3.fa");

    if (CO.verbose)
    {
      std::cout << "Adding utr    " << base_path << std::endl;
    }

    // graph = createGraph(base_path.c_str(), vertex_labels, ids, begin_vertex);
    return;
  }
  else if (p5_region)
  {
    std::string utr = "5p";
    utr.append(".fa");
    base_path.append(utr);

    if (CO.verbose)
    {
      std::cout << "Adding utr    " << base_path << std::endl;
    }
  }
  else if (intron_region)
  {
    std::string intron = "i";
    std::string feature_number_str = std::to_string(feature_number);
    intron.append(feature_number_str);
    intron.append(".fa");
    base_path.append(intron);

    if (CO.verbose)
    {
      std::cout << "Adding intron " << base_path << std::endl;
    }
  }
  else
  {
    // Exon region
    std::string exon = "e";
    std::string feature_number_str = std::to_string(feature_number);
    exon.append(feature_number_str);
    exon.append(".fa");
    base_path.append(exon);
  }

  free_nodes.insert(begin_vertex);

  if (add_bitstrings)
  {
    // extendGraph(graph, base_path.c_str(), vertex_labels, edge_ids, new_begin_vertex, begin_vertex);
  }
  else
  {
    // extendGraph(graph, base_path.c_str(), vertex_labels, new_begin_vertex, begin_vertex);
  }
}

bool
Constructor::read_reference_genome(const char * fasta_filename)
{
  if (!seqan::open(fasta_index, fasta_filename))
  {
    if (!seqan::build(fasta_index, fasta_filename))
    {
      std::cerr << "ERROR: Index could not be loaded or built." << std::endl;
      return false;
    }

    if (!seqan::save(fasta_index))
    {
      std::cerr << "ERROR: Index could not be written do disk." << std::endl;
      return false;
    }
  }

  return true;
}

bool
Constructor::extract_reference_sequence(const char * region)
{
  return true;
}

unsigned
Constructor::get_fasta_index_id(const char * id)
{
  SEQAN_CHECK(numSeqs(fasta_index) > 0, "No sequences found in FASTA index.");
  unsigned fasta_index_id = 0;
  int ret = seqan::getIdByName(fasta_index_id, fasta_index, id);
  SEQAN_CHECK(ret, "Could not find the FASTA id = %s", id);
  return fasta_index_id;
}

void
Constructor::open_vcf(const char * vcf_filename)
{
  seqan::open(vcf_file, vcf_filename);
  seqan::VcfHeader header;
  seqan::readHeader(header, vcf_file);
}

int
Constructor::read_vcf_record()
{
  if (!seqan::atEnd(vcf_file))
  {
    seqan::readRecord(vcf_record, vcf_file);
    return 0;
  }

  return 1;
}

void
Constructor::open_tabix(const char * tabix_filename)
{
  seqan::open(tabix_file, tabix_filename);
  seqan::String<char> header;
  seqan::getHeader(header, tabix_file);
  // std::cout << "header = " << header << std::endl;
}

bool
Constructor::read_tabix_record()
{
  return seqan::readRecord(vcf_record, tabix_file);
}


} // namespace gyper

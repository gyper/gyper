#include "partial_order_graph.hpp"


Gyper::Gyper ()
{
	TGraph graph();
}

Gyper::Gyper (Options & CO)
{
  TGraph graph();
  Gyper::Gyper::CO = CO;
}

void
Gyper::create_HLA_graph()
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
  topologicalSort(order, graph);
}

void
Gyper::index()
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
Gyper::get_number_of_exons()
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
Gyper::get_HLA_base_path()
{
  std::stringstream base_path;
  base_path << gyper_SOURCE_DIRECTORY << "/data/haplotypes/hla/references/" << CO.gene << "/";
  return base_path.str();
}

void
Gyper::add_FASTA_region(bool add_bitstrings, int feature_number, bool intron_region, bool p3_region, bool p5_region)
{
  std::string base_path = get_HLA_base_path();

  if (p3_region)
  {
    base_path.append("p3.fa");

    if (CO.verbose)
    {
      std::cout << "Adding utr    " << base_path << std::endl;
    }

    graph = createGraph(base_path.c_str(), vertex_vector, ids, begin_vertex);
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
    extendGraph(graph, base_path.c_str(), vertex_vector, edge_ids, new_begin_vertex, begin_vertex);
  }
  else
  {
    extendGraph(graph, base_path.c_str(), vertex_vector, new_begin_vertex, begin_vertex);
  }
}

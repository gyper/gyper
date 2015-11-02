#include "partial_order_graph.hpp"

unsigned get_number_of_exons(seqan::String<char> gene)
{
  // Set the number of exons
  if (gene == "HLAA" || gene == "HLAC")
  {
    return 8;
  }
  else if (gene == "HLAB")
  {
    return 7;
  }
  else if (gene == "DQA1")
  {
    return 4;
  }
  else if (gene == "DQB1" || gene == "DRB1")
  {
    return 6;
  }
  std::cerr << "Unsupported HLA gene: " << gene << std::endl;
  return 0;
}

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
  unsigned number_of_exons = get_number_of_exons(CO.gene);

  std::stringstream base_path;
  base_path << gyper_SOURCE_DIRECTORY << "/data/haplotypes/hla/references/" << CO.gene << "/";
  TVertexDescriptor begin_vertex;

  {
    std::string p3_string = base_path.str();
    p3_string.append("p3.fa");

    if (CO.verbose)
    {
      std::cout << "Adding utr    " << p3_string << std::endl;
    }

    const char* alignment_file_p3 = p3_string.c_str();
    graph = createGraph(alignment_file_p3, vertex_vector, ids, begin_vertex);
  }
  
  TVertexDescriptor new_begin_vertex;
  std::string tmp_string;
  std::string extension = ".fa";

  // First exon
  {
    tmp_string = base_path.str();
    std::string exon = "e";
    std::string feature_number = std::to_string(number_of_exons);
    exon.append(feature_number);
    exon.append(extension);
    tmp_string.append(exon);

    // if (CO.verbose)
    // {
    //   std::cout << "Adding exon   " << tmp_string << std::endl;
    // }

    const char* alignment_file = tmp_string.c_str();
    free_nodes.insert(begin_vertex);

    if (CO.verbose)
    {
      std::cout << "Adding exon   " << tmp_string << " as intron" << std::endl;
    }

    extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);

    --number_of_exons;
  }

  while (number_of_exons >= 1)
  {
    // Intron
    {
      tmp_string = base_path.str();
      std::string intron = "i";
      std::string feature_number = std::to_string(number_of_exons);
      intron.append(feature_number);
      intron.append(extension);
      tmp_string.append(intron);

      if (CO.verbose)
      {
        std::cout << "Adding intron " << tmp_string << std::endl;
      }

      const char* alignment_file = tmp_string.c_str();
      free_nodes.insert(begin_vertex);
      extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);
    }

    // Exon
    {
      tmp_string = base_path.str();
      std::string exon = "e";
      std::string feature_number = std::to_string(number_of_exons);
      exon.append(feature_number);
      exon.append(extension);
      tmp_string.append(exon);

      // if (CO.verbose)
      // {
      //   std::cout << "Adding exon   " << tmp_string << std::endl;
      // }

      const char* alignment_file = tmp_string.c_str();
      free_nodes.insert(begin_vertex);
      
      if (CO.exon_2_and_3 && (number_of_exons == 2 || number_of_exons == 3 || number_of_exons == 4))
      {
        if (CO.verbose)
        {
          std::cout << "Adding exon   " << tmp_string << std::endl;
        }

        extendGraph(graph, alignment_file, vertex_vector, edge_ids, new_begin_vertex, begin_vertex);
      }
      else
      {
        if (CO.verbose)
        {
          std::cout << "Adding exon   " << tmp_string << " as intron" << std::endl;
        }

        extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);
      }
    }

    --number_of_exons;
  }

  {
    // Final UTR
    tmp_string = base_path.str();
    std::string utr = "5p";
    utr.append(extension);
    tmp_string.append(utr);

    if (CO.verbose)
    {
      std::cout << "Adding utr    " << tmp_string << std::endl;
    }

    const char* alignment_file = tmp_string.c_str();
    free_nodes.insert(begin_vertex);
    extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);
  }

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

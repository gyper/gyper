#include "graph.hpp"


void
align_sequence (DnaString & my_sequence,
                boost::dynamic_bitset<> & qual,
                TGraph const & graph,
                std::vector<VertexLabels> & vertex_vector,
                String<TVertexDescriptor> & order,
                std::vector<ExactBacktracker> & backtracker,
                std::vector<ExactBacktracker> & reverse_backtracker,
                boost::unordered_set<TVertexDescriptor> const & free_nodes,
                std::vector<TVertexDescriptor> & matching_vertices,
                std::vector<TVertexDescriptor> & reverse_matching_vertices
               )
{
  initializeExactScoreMatrixAndBacktracker(length(my_sequence), length(order), backtracker);
  reverse_backtracker = backtracker;

  alignToGraphExact (my_sequence,
                     order,
                     graph,
                     matching_vertices,
                     vertex_vector,
                     backtracker,
                     free_nodes,
                     qual
                    );

  reverseComplement(my_sequence);

  boost::dynamic_bitset<> qual_reversed(qual.size());
  std::size_t qual_size = qual.size();
  for (unsigned pos = 0 ; pos < qual_size ; ++pos)
  {
    if (qual.test(pos))
    {
      qual_reversed[qual_size-pos-1] = 1;
    }
  }
  
  alignToGraphExact (my_sequence,
                     order,
                     graph,
                     reverse_matching_vertices,
                     vertex_vector,
                     reverse_backtracker,
                     free_nodes,
                     qual
                    );

  reverseComplement(my_sequence);
}


boost::dynamic_bitset<>
align_sequence_kmer (String<Dna> & my_sequence,
                     String<char> & qual,
                     unsigned const & id_numbers,
                     TKmerMap & kmer_map,
                     std::vector<VertexLabels> & vertex_vector,
                     int const & kmer_size,
                     int const & min_kmers
                    )
{
  unsigned best_kmer_index = find_best_kmer(qual, kmer_size);
  boost::dynamic_bitset<> matched_ids =
   align_kmer_to_graph (my_sequence,
                        id_numbers,
                        kmer_map,
                        vertex_vector,
                        best_kmer_index,
                        kmer_size,
                        min_kmers
                       );

  if (matched_ids.find_first() != matched_ids.npos)
  {
    return matched_ids;
  }

  reverseComplement(my_sequence);
  matched_ids =
   align_kmer_to_graph (my_sequence,
                        id_numbers,
                        kmer_map,
                        vertex_vector,
                        best_kmer_index,
                        kmer_size,
                        min_kmers
                       );

  reverseComplement(my_sequence);
  return matched_ids;
}


CharString myExtractTagValue(String<char> &tags)
{
  BamTagsDict tagsDict(tags);
  unsigned tagIdx = 0;
  if (!findTagKey(tagIdx, tagsDict, "RG"))
  {
    return "";
  }

  CharString read_group;

  if (!extractTagValue(read_group, tagsDict, tagIdx))
  {
    std::cerr << "Not a valid string at pos " << tagIdx << std::endl;
    return "";
  }

  return read_group;
}


void
createGenericGraph(callOptions & CO,
                   TGraph & graph,
                   std::vector<VertexLabels> & vertex_vector,
                   std::vector<std::string> & ids,
                   boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                   boost::unordered_set<TVertexDescriptor> & free_nodes,
                   String<TVertexDescriptor> & order
                  )
{
  int number_of_exons = CO.number_of_exons;

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

    if (CO.verbose)
    {
      std::cout << "Adding exon   " << tmp_string << std::endl;
    }

    const char* alignment_file = tmp_string.c_str();
    free_nodes.insert(begin_vertex);
    extendGraph(graph, alignment_file, vertex_vector, edge_ids, new_begin_vertex, begin_vertex);

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

      if (CO.verbose)
      {
        std::cout << "Adding exon   " << tmp_string << std::endl;
      }

      const char* alignment_file = tmp_string.c_str();
      free_nodes.insert(begin_vertex);
      extendGraph(graph, alignment_file, vertex_vector, edge_ids, new_begin_vertex, begin_vertex);
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
create_exon_2_and_3_graph(callOptions & CO,
                          TGraph & graph,
                          std::vector<VertexLabels> & vertex_vector,
                          std::vector<std::string> & ids,
                          boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                          boost::unordered_set<TVertexDescriptor> & free_nodes,
                          String<TVertexDescriptor> & order
                         )
{
  int number_of_exons = CO.number_of_exons;

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
    if (number_of_exons == 2 || number_of_exons == 3)
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
      if (number_of_exons == 2 || number_of_exons == 3)
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

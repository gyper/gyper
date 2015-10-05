#include "graph_kmerify.hpp"


boost::unordered_map<String<Dna>, KmerLabels >
kmerifyGraph(String<TVertexDescriptor const> const & order,
             TGraph const & graph,
             std::vector<VertexLabels> & vertex_vector,
             boost::unordered_set<TVertexDescriptor> const & free_nodes
            )
{
  boost::unordered_map<String<Dna>, KmerLabels > kmer_map;
  std::vector< std::vector<KmerLabelsWithKmer> > kmer_tracker;
  std::vector<KmerLabelsWithKmer> kmer_labels;
  KmerLabelsWithKmer kmer_label;
  kmer_labels.push_back(kmer_label);
  kmer_tracker.resize(length(order), kmer_labels);

  Iterator<String<TVertexDescriptor const> const>::Type it = begin(order);

  std::cout << "Processing first source." << std::endl;
  {
    // First source vertex is special
    TVertexDescriptor const & first_source_vertex = *it;
    std::cout << "first_source_vertex, length(order) = " << first_source_vertex << " " << length(order) << std::endl;
    KmerLabelsWithKmer source_label = kmer_tracker[first_source_vertex][0];
    std::cout << "Iterating outer edges." << std::endl;

    for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, first_source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
    {
      TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);
      VertexLabels target_vertex_label = vertex_vector[target_vertex];
      KmerLabelsWithKmer target_label = kmer_tracker[target_vertex][0];
      target_label.kmer = target_vertex_label.dna;
      target_label.nodes.push_back(target_vertex);
      target_label.edges.flip();
      std::cout << "Target label " << target_vertex << " is now: " << target_label << std::endl;
    }
  }

  for ( ; it != end(order) ; ++it)
  {
    TVertexDescriptor const & source_vertex = *it;

    for (unsigned index_of_source_labels = 0 ; index_of_source_labels < kmer_tracker[source_vertex].size() ; ++index_of_source_labels)
    {
      KmerLabelsWithKmer const & source_label = kmer_tracker[source_vertex][index_of_source_labels];

      for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
      {
        if (length(source_label.kmer) == K_SIZE)
        {
          // If the current kmer is of length K, then we must first remove the first base
          // if (kmer_map.find(source_label.kmer) != kmer_map.end())
          // {
          //   std::cout << "Added kmer " << source_label.kmer << std::endl;
          //   KmerLabels new_kmer_label(source_vertex);
          //   kmer_map[source_label.kmer] = new_kmer_label;
          // }
          // source_label.kmer = 
        }

        TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);

        if (free_nodes.find(target_vertex) != free_nodes.end())
        {
          // Target is a "free node"

        }
        else
        {
          // Target is a typical node

        }
      }
    }
  }

  return kmer_map;
}

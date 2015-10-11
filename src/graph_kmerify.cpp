#include "graph_kmerify.hpp"


void
checkKmers(DnaString const & kmer,
           TVertexDescriptor const & source_vertex,
           TGraph const & graph,
           std::vector<VertexLabels> & vertex_vector,
           boost::unordered_set<TVertexDescriptor> const & free_nodes,
           TKmerMap & kmer_map)
{
  if (length(kmer) == K_SIZE)
  {
    if (kmer_map.count(kmer) == 0)
    {
      std::vector<TVertexDescriptor> new_vector(1, source_vertex);
      kmer_map[kmer] = new_vector;
      // std::cout << "Added " << kmer << std::endl;
    }
    else
    {
      kmer_map[kmer].push_back(source_vertex);
    }

    return;
  }

  for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
  {
    DnaString new_kmer(kmer);
    TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);

    if (free_nodes.count(target_vertex) == 0)
    {
      seqan::appendValue(new_kmer, vertex_vector[target_vertex].dna);
    }

    checkKmers(new_kmer, target_vertex, graph, vertex_vector, free_nodes, kmer_map);
  }
}


TKmerMap
kmerifyGraph(String<TVertexDescriptor const> const & order,
             TGraph const & graph,
             std::vector<VertexLabels> & vertex_vector,
             boost::unordered_set<TVertexDescriptor> const & free_nodes
            )
{
  TKmerMap kmer_map;

  for (Iterator<String<TVertexDescriptor const> const>::Type it = begin(order) ; it != end(order) ; ++it)
  {
    TVertexDescriptor const & source_vertex = *it;
    checkKmers("", source_vertex, graph, vertex_vector, free_nodes, kmer_map);
    std::cout << "kmer_map.size() = " << kmer_map.size() << std::endl;
  }

  return kmer_map;
}

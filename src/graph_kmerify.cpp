#include "graph_kmerify.hpp"


void
checkKmers(DnaString const & kmer,
           TVertexDescriptor const & starting_vertex,
           TVertexDescriptor const & source_vertex,
           TGraph const & graph,
           std::vector<VertexLabels> & vertex_vector,
           boost::unordered_set<TVertexDescriptor> const & free_nodes,
           TKmerMap & kmer_map,
           int const & kmer_size
          )
{
  if (length(kmer) == static_cast<std::size_t>(kmer_size))
  {
    if (kmer_map.count(kmer) == 0)
    {
      std::vector<TVertexDescriptor> new_vector(1, starting_vertex);
      kmer_map[kmer] = new_vector;
      // std::cout << "Added " << kmer << std::endl;
    }
    else
    {
      std::vector<TVertexDescriptor>::iterator vertex_it = find (kmer_map[kmer].begin(), kmer_map[kmer].end(), starting_vertex);

      if (vertex_it == kmer_map[kmer].end())
      {
        kmer_map[kmer].push_back(starting_vertex);
      }
    }

    return;
  }

  for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
  {
    DnaString new_kmer(kmer);
    TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);

    // std::cout << source_vertex << " -> " << target_vertex << std::endl;

    if (free_nodes.count(target_vertex) == 0)
    {
      seqan::appendValue(new_kmer, vertex_vector[target_vertex].dna);
    }

    checkKmers(new_kmer, starting_vertex, target_vertex, graph, vertex_vector, free_nodes, kmer_map, kmer_size);
  }
}


TKmerMap
kmerifyGraph(String<TVertexDescriptor const> const & order,
             TGraph const & graph,
             std::vector<VertexLabels> & vertex_vector,
             boost::unordered_set<TVertexDescriptor> const & free_nodes,
             int const & kmer_size
            )
{
  TKmerMap kmer_map;

  for (Iterator<String<TVertexDescriptor const> const>::Type it = begin(order) ; it != end(order) ; ++it)
  {
    TVertexDescriptor const & source_vertex = *it;

    if (free_nodes.count(source_vertex) == 0)
    {
      checkKmers(vertex_vector[source_vertex].dna, source_vertex, source_vertex, graph, vertex_vector, free_nodes, kmer_map, kmer_size);
      // std::cout << "source_vertex = " << source_vertex << " kmer_map.size() = " << kmer_map.size() << " kmer_size = " << kmer_size;
      // std::cout << " vertex_vector[source_vertex].level = " << vertex_vector[source_vertex].level << std::endl;
    }
  }

  return kmer_map;
}

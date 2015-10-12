#include "graph_kmerify.hpp"


void
checkKmers(DnaString const & kmer,
           TVertexDescriptor const & starting_vertex,
           TVertexDescriptor const & source_vertex,
           TGraph const & graph,
           std::vector<VertexLabels> & vertex_vector,
           boost::unordered_set<TVertexDescriptor> const & free_nodes,
           boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
           boost::dynamic_bitset<> const & id_bits,
           TKmerMap & kmer_map,
           std::size_t const & kmer_size
          )
{
  if (length(kmer) == kmer_size)
  {
    KmerLabels new_kmer_label = {
      starting_vertex,
      source_vertex,
      id_bits
    };

    if (kmer_map.count(kmer) == 0)
    {
      std::vector<KmerLabels> new_vector(1, new_kmer_label);
      kmer_map[kmer] = new_vector;
    }
    else
    {
      // for (std::vector<KmerLabels>::iterator vertex_it = kmer_map[kmer].begin() ; vertex_it != kmer_map[kmer].end() ; ++vertex_it)
      // {
      //   if (vertex_it->start_vertex == starting_vertex)
      //   {
      //     return;
      //   }
      // }

      kmer_map[kmer].push_back(new_kmer_label);
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
      std::pair<TVertexDescriptor, TVertexDescriptor> edge_pair(source_vertex, target_vertex);
      
      boost::dynamic_bitset<> new_id_bits(id_bits);

      if (edge_ids.count(edge_pair) == 1)
      {
        // std::cout << new_id_bits << " to ";
        new_id_bits = id_bits & edge_ids[edge_pair];
        // std::cout << new_id_bits << " (" << edge_ids[edge_pair] << ")" << std::endl;
      }

      checkKmers(new_kmer, starting_vertex, target_vertex, graph, vertex_vector, free_nodes, edge_ids, new_id_bits, kmer_map, kmer_size);
    }
    else
    {
      checkKmers(new_kmer, starting_vertex, target_vertex, graph, vertex_vector, free_nodes, edge_ids, id_bits, kmer_map, kmer_size);
    }
  }
}


TKmerMap
kmerifyGraph(String<TVertexDescriptor const> const & order,
             TGraph const & graph,
             std::vector<VertexLabels> & vertex_vector,
             boost::unordered_set<TVertexDescriptor> const & free_nodes,
             boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
             int const & kmer_size
            )
{
  TKmerMap kmer_map;

  for (Iterator<String<TVertexDescriptor const> const>::Type it = begin(order) ; it != end(order) ; ++it)
  {
    TVertexDescriptor const & source_vertex = *it;

    if (free_nodes.count(source_vertex) == 0)
    {
      boost::dynamic_bitset<> id_bits(edge_ids.begin()->second.size());
      id_bits.flip();
      checkKmers(vertex_vector[source_vertex].dna, source_vertex, source_vertex, graph, vertex_vector, free_nodes, edge_ids, id_bits, kmer_map, static_cast<std::size_t>(kmer_size));
      // std::cout << "source_vertex = " << source_vertex << " kmer_map.size() = " << kmer_map.size() << " kmer_size = " << kmer_size;
      // std::cout << " vertex_vector[source_vertex].level = " << vertex_vector[source_vertex].level << std::endl;
    }
  }

  return kmer_map;
}

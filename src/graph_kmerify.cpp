#include "graph_kmerify.hpp"

#include <cmath>

using namespace seqan;

void
check_kmers_simple(DnaString const & kmer,
                   TGraph const & graph,
                   TVertexDescriptor const & source_vertex,
                   std::vector<VertexLabels> & vertex_vector,
                   boost::unordered_set<TVertexDescriptor> const & free_nodes,
                   boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                   boost::dynamic_bitset<> const & id_bits,
                   TKmerMapSimple & kmer_map,
                   std::size_t const & kmer_size
                  )
{
  if (id_bits.none())
    return;

  if (length(kmer) == kmer_size)
  {
    if (id_bits.all())
      return;

    if (kmer_map.count(kmer) == 0)
    {
      kmer_map[kmer] = id_bits;
    }
    else
    {
      kmer_map[kmer] |= id_bits;
    }

    return;
  }

  for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
  {
    DnaString new_kmer(kmer);
    TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);

    // std::cout << source_vertex << " -> " << target_vertex << std::endl;

    boost::dynamic_bitset<> new_id_bits(id_bits);

    if (free_nodes.count(target_vertex) == 0)
    {
      seqan::appendValue(new_kmer, vertex_vector[target_vertex].dna);
      std::pair<TVertexDescriptor, TVertexDescriptor> edge_pair(source_vertex, target_vertex);
      
      if (edge_ids.count(edge_pair) == 1)
      {
        new_id_bits = id_bits & edge_ids[edge_pair];
      }
    }

    check_kmers_simple(new_kmer, graph, target_vertex, vertex_vector, free_nodes, edge_ids, new_id_bits, kmer_map, kmer_size);
  }
}


TKmerMapSimple
kmerify_graph_simple(String<TVertexDescriptor const> const & order,
                     TGraph const & graph,
                     std::vector<VertexLabels> & vertex_vector,
                     boost::unordered_set<TVertexDescriptor> const & free_nodes,
                     boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                     int const & kmer_size
                    )
{
  TKmerMapSimple kmer_map;

  for (Iterator<String<TVertexDescriptor const> const>::Type it = begin(order) ; it != end(order) ; ++it)
  {
    TVertexDescriptor const & source_vertex = *it;

    if (free_nodes.count(source_vertex) == 0)
    {
      boost::dynamic_bitset<> id_bits(edge_ids.begin()->second.size());
      id_bits.flip();
      check_kmers_simple(vertex_vector[source_vertex].dna, graph, source_vertex, vertex_vector, free_nodes, edge_ids, id_bits, kmer_map, static_cast<std::size_t>(kmer_size));
    }
  }

  return kmer_map;
}



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
  if (id_bits.none())
    return;

  if (length(kmer) == kmer_size)
  {
    KmerLabels new_kmer_label =
    {
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
      kmer_map[kmer].push_back(new_kmer_label);
    }

    return;
  }

  for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
  {
    DnaString new_kmer(kmer);
    TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);
    boost::dynamic_bitset<> new_id_bits(id_bits);

    if (free_nodes.count(target_vertex) == 0)
    {
      seqan::appendValue(new_kmer, vertex_vector[target_vertex].dna);
      std::pair<TVertexDescriptor, TVertexDescriptor> edge_pair(source_vertex, target_vertex);
      
      if (edge_ids.count(edge_pair) == 1)
      {
        new_id_bits = id_bits & edge_ids[edge_pair];
      }
    }

    checkKmers(new_kmer, starting_vertex, target_vertex, graph, vertex_vector, free_nodes, edge_ids, new_id_bits, kmer_map, kmer_size);
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


unsigned
find_best_kmer(String<char> qual,
               unsigned const & k
              )
{
  unsigned best_index = 0;
  double P_error = 1;

  auto end_it = begin(qual);
  auto first_end = begin(qual)+k;
  
  for ( ; end_it != first_end ; ++end_it)
  {
    P_error *= 1.0 - std::pow(10.0, static_cast<double>(*end_it - 33) / -10.0);
  }

  double min_P_error = P_error;
  unsigned index = 1;

  for (auto start_it = begin(qual) ; end_it != end(qual) ; ++start_it, ++end_it, ++index)
  {
    // std::cout << P_error << std::endl;
    P_error *= 1.0 - std::pow(10.0, static_cast<double>(*end_it - 33) / -10.0);
    P_error /= 1.0 - std::pow(10.0, static_cast<double>(*start_it - 33) / -10.0);
    
    if (P_error > min_P_error)
    {
      best_index = index;
      min_P_error = P_error;
    }
  }

  // std::cout << P_error << std::endl;
  return best_index;
}

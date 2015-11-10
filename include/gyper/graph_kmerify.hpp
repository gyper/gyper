#ifndef __GRAPH_KMERIFY_HPP_INCLUDED__
#define __GRAPH_KMERIFY_HPP_INCLUDED__
#include <string>
#include <sstream>
#include <iostream>

#include <gyper/graph.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>


TKmerMapSimple
kmerify_graph_simple(seqan::String<TVertexDescriptor const> const & order,
                     TGraph const & graph,
                     std::vector<VertexLabels> & vertex_vector,
                     boost::unordered_set<TVertexDescriptor> const & free_nodes,
                     boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                     int const & kmer_size
                    );


TKmerMap
kmerifyGraph(seqan::String<TVertexDescriptor const> const & order,
             TGraph const & graph,
             std::vector<VertexLabels> & vertex_vector,
             boost::unordered_set<TVertexDescriptor> const & free_nodes,
             boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
             int const & kmer_size
            );


unsigned
find_best_kmer(seqan::String<char> qual,
               unsigned const & k
              );

#endif

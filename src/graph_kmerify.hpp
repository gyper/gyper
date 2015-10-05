#ifndef __GRAPH_KMERIFY_HPP_INCLUDED__
#define __GRAPH_KMERIFY_HPP_INCLUDED__
#define K_SIZE 15
#include <string>
#include <sstream>
#include <iostream>

#include "graph_builder.hpp"
#include "graph_align.hpp"

struct KmerLabels {
  int start_node;
  boost::dynamic_bitset<> edges;

  KmerLabels(int starting_node) : start_node(starting_node),  edges(K_SIZE) {}

};

static inline std::ostream &
operator<< (std::ostream& os, KmerLabels const & rhs)
{
  os << rhs.start_node << " " << rhs.edges;
  return os;
}

struct KmerLabelsWithKmer {
  String<Dna> kmer;
  std::list<int> nodes;
  boost::dynamic_bitset<> edges;

  KmerLabelsWithKmer() : kmer(""), nodes(), edges(K_SIZE) {}
};

static inline std::ostream &
operator<< (std::ostream& os, KmerLabelsWithKmer const & rhs)
{
  os << rhs.kmer << " " << rhs.nodes.size() << " " << rhs.edges;
  return os;
}

boost::unordered_map<String<Dna>, KmerLabels >
kmerifyGraph(String<TVertexDescriptor const> const & order,
             TGraph const & graph,
             std::vector<VertexLabels> & vertex_vector,
             boost::unordered_set<TVertexDescriptor> const & free_nodes
            );

#endif

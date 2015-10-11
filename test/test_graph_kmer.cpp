#include <stdio.h>
#include <climits>

#include "../src/graph_kmerify.hpp"
#include "constants.hpp"
#include "catch.hpp"

// Use 4-mers when testing
#undef K_SIZE
#define K_SIZE 4 


using namespace seqan;

TEST_CASE ("Kmerify should work")
{
  TGraph g;
  CharString sequence = "TAGAT";
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  std::vector<VertexLabels> vertex_vector;
  TVertexDescriptor begin_vertex = 0;

  addInitialAndEndVertex(g, vertex_label_map, vertex_vector, begin_vertex);
  addExonToGraph(g, 0, 1, sequence, vertex_label_map, vertex_vector, edge_ids);
}


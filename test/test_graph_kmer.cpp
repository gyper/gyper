#include <stdio.h>
#include <climits>

#define TEST_GYPER

#include "../src/graph_kmerify.hpp"
#include "constants.hpp"
#include "catch.hpp"

using namespace seqan;

TEST_CASE ("kmerifyGraph should create a non-empty kmer map")
{
  std::ostringstream ss;
  ss << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case5_e3.fa";
  std::string tmp_str = ss.str();
  const char* alignment_file_e3 = tmp_str.c_str();

  std::vector<std::string> ids;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;
  boost::unordered_set<TVertexDescriptor> free_nodes;
  free_nodes.insert(1);
  free_nodes.insert(0);

  TVertexDescriptor begin_vertex;
  TGraph graph = createGraph(alignment_file_e3, vertex_vector, edge_ids, ids, begin_vertex);

  // Use 4-mers when testing
  int kmer_size = 4;

  SECTION("We should only find one kmer for test_case5, exon 3")
  {
    String<TVertexDescriptor> order;
    topologicalSort(order, graph);

    TKmerMap kmer_map = kmerifyGraph(order, graph, vertex_vector, free_nodes, kmer_size);
    REQUIRE(kmer_map.size() == 1);
    DnaString kmer1 = "TGCC";
    REQUIRE(kmer_map.count(kmer1) == 1);
    REQUIRE(kmer_map[kmer1].size() == 1);
    REQUIRE(kmer_map[kmer1][0] == 2);
  }

  std::ostringstream ss2;
  ss2 << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case5_i2.fa";
  std::string tmp_str2 = ss2.str();
  const char* alignment_file_i2 = tmp_str2.c_str();
  TVertexDescriptor new_begin_vertex;
  free_nodes.insert(begin_vertex);
  extendGraph(graph, alignment_file_i2, vertex_vector, new_begin_vertex, begin_vertex);
  std::vector<ExactBacktracker> backtracker;

  
  SECTION ("Add intron 2")
  {
    String<TVertexDescriptor> order;
    topologicalSort(order, graph);
    TKmerMap kmer_map = kmerifyGraph(order, graph, vertex_vector, free_nodes, kmer_size);
    // REQUIRE(kmer_map.size() == 4);

    DnaString kmer1 = "TTTT";
    REQUIRE(kmer_map.count(kmer1) == 1);
    REQUIRE(kmer_map[kmer1].size() == 2);
    REQUIRE(kmer_map[kmer1][0] == 7);
    REQUIRE(kmer_map[kmer1][1] == 8);

    DnaString kmer2 = "TTTG";
    REQUIRE(kmer_map.count(kmer2) == 1);
    REQUIRE(kmer_map[kmer2].size() == 1);
    REQUIRE(kmer_map[kmer2][0] == 9);

    DnaString kmer3 = "TTGC";
    REQUIRE(kmer_map.count(kmer3) == 1);
    REQUIRE(kmer_map[kmer3].size() == 1);
    REQUIRE(kmer_map[kmer3][0] == 10);

    DnaString kmer4 = "TGCC";
    REQUIRE(kmer_map.count(kmer4) == 1);
    REQUIRE(kmer_map[kmer4].size() == 1);
    REQUIRE(kmer_map[kmer4][0] == 2);

    /*
    SECTION ("sequence is GATA")
    {
      DnaString sequence = "GATA";

      initializeExactScoreMatrixAndBacktracker(length(sequence), length(order), backtracker);

      std::vector<TVertexDescriptor> matching_vertices;
      alignToGraphExact (sequence,
                       order,
                       graph,
                       matching_vertices,
                       vertex_vector,
                       backtracker,
                       free_nodes
                      );

      REQUIRE (matching_vertices.size() == 1);
      REQUIRE (matching_vertices[0] == 5);

      boost::dynamic_bitset<> perfect_match =
        backTrackAndCount(ids, backtracker, matching_vertices[0], edge_ids);

      boost::dynamic_bitset<> expected(2, 1ul); // Binary: 01
      REQUIRE (expected == perfect_match);
    }

    SECTION ("sequence is AGATA")
    {
      DnaString sequence = "AGATA";

      initializeExactScoreMatrixAndBacktracker(length(sequence), length(order), backtracker);

      std::vector<TVertexDescriptor> matching_vertices;
      alignToGraphExact (sequence,
                       order,
                       graph,
                       matching_vertices,
                       vertex_vector,
                       backtracker,
                       free_nodes
                      );

      REQUIRE (matching_vertices.size() == 1);
      REQUIRE (matching_vertices[0] == 5);

      boost::dynamic_bitset<> perfect_match =
        backTrackAndCount(ids, backtracker, matching_vertices[0], edge_ids);

      boost::dynamic_bitset<> expected(2, 1ul); // Binary: 00001
      REQUIRE (expected == perfect_match);
    }
    */
  }
}

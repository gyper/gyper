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


  SECTION("We should only find one 4 kmer for test_case5, exon 3")
  {
    int kmer_size = 4;
    String<TVertexDescriptor> order;
    topologicalSort(order, graph);

    TKmerMap kmer_map = kmerifyGraph(order, graph, vertex_vector, free_nodes, edge_ids, kmer_size);
    REQUIRE(kmer_map.size() == 2);

    DnaString kmer1 = "TGCC";
    REQUIRE(kmer_map.count(kmer1) == 1);
    REQUIRE(kmer_map[kmer1].size() == 1);
    REQUIRE(kmer_map[kmer1][0].start_vertex == 2);
    REQUIRE(kmer_map[kmer1][0].end_vertex == 5);

    DnaString kmer2 = "TGCA";
    REQUIRE(kmer_map.count(kmer2) == 1);
    REQUIRE(kmer_map[kmer2].size() == 1);
    REQUIRE(kmer_map[kmer2][0].start_vertex == 2);
    REQUIRE(kmer_map[kmer2][0].end_vertex == 6);
  }

  std::ostringstream ss2;
  ss2 << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case5_i2.fa";
  std::string tmp_str2 = ss2.str();
  const char* alignment_file_i2 = tmp_str2.c_str();
  TVertexDescriptor new_begin_vertex;
  free_nodes.insert(begin_vertex);
  extendGraph(graph, alignment_file_i2, vertex_vector, new_begin_vertex, begin_vertex);
  free_nodes.insert(begin_vertex);

  std::size_t num_ids = ids.size();
  

  SECTION ("Add intron 2 and align using 4 mers")
  {
    // Use 4-mers when testing
    int kmer_size = 4;
    String<TVertexDescriptor> order;
    topologicalSort(order, graph);
    TKmerMap kmer_map = kmerifyGraph(order, graph, vertex_vector, free_nodes, edge_ids, kmer_size);

    boost::dynamic_bitset<> expected_bitset_11(num_ids, 3ul);
    boost::dynamic_bitset<> expected_bitset_10(num_ids, 2ul);
    boost::dynamic_bitset<> expected_bitset_01(num_ids, 1ul);
    boost::dynamic_bitset<> expected_bitset_00(num_ids, 0ul);

    REQUIRE(kmer_map.size() == 5);

    DnaString kmer1 = "TTTT";
    REQUIRE(kmer_map.count(kmer1) == 1);
    REQUIRE(kmer_map[kmer1].size() == 2);
    REQUIRE(kmer_map[kmer1][0].start_vertex == 8);
    REQUIRE(kmer_map[kmer1][0].end_vertex == 11);
    REQUIRE(kmer_map[kmer1][0].id_bits == expected_bitset_11);
    REQUIRE(kmer_map[kmer1][1].start_vertex == 9);
    REQUIRE(kmer_map[kmer1][1].end_vertex == 2);
    
    REQUIRE(kmer_map[kmer1][1].id_bits == expected_bitset_11);

    DnaString kmer2 = "TTTG";
    REQUIRE(kmer_map.count(kmer2) == 1);
    REQUIRE(kmer_map[kmer2].size() == 1);
    REQUIRE(kmer_map[kmer2][0].start_vertex == 10);
    REQUIRE(kmer_map[kmer2][0].end_vertex == 3);
    REQUIRE(kmer_map[kmer2][0].id_bits == expected_bitset_11);

    DnaString kmer3 = "TTGC";
    REQUIRE(kmer_map.count(kmer3) == 1);
    REQUIRE(kmer_map[kmer3].size() == 2);
    REQUIRE(kmer_map[kmer3][0].start_vertex == 11);
    REQUIRE(kmer_map[kmer3][0].end_vertex == 5);
    REQUIRE(kmer_map[kmer3][0].id_bits == expected_bitset_10);
    REQUIRE(kmer_map[kmer3][1].start_vertex == 11);
    REQUIRE(kmer_map[kmer3][1].end_vertex == 4);
    REQUIRE(kmer_map[kmer3][1].id_bits == expected_bitset_01);

    DnaString kmer4 = "TGCC";
    REQUIRE(kmer_map.count(kmer4) == 1);
    REQUIRE(kmer_map[kmer4].size() == 1);
    REQUIRE(kmer_map[kmer4][0].start_vertex == 2);
    REQUIRE(kmer_map[kmer4][0].end_vertex == 5);
    REQUIRE(kmer_map[kmer4][0].id_bits == expected_bitset_01);

    DnaString kmer5 = "TGCA";
    REQUIRE(kmer_map.count(kmer5) == 1);
    REQUIRE(kmer_map[kmer5].size() == 1);
    REQUIRE(kmer_map[kmer5][0].start_vertex == 2);
    REQUIRE(kmer_map[kmer5][0].end_vertex == 6);
    REQUIRE(kmer_map[kmer5][0].id_bits == expected_bitset_10);

    DnaString kmer6 = "GCCA"; // Not possible! Since the bitstring will be 00
    REQUIRE(kmer_map.count(kmer6) == 0);

    SECTION ("sequence is TGCC")
    {
      DnaString sequence = "TGCC";
      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_01 == actual_bitset);
    }

    SECTION ("sequence is TTGC")
    {
      DnaString sequence = "TTGC";
      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_11 == actual_bitset);
    }

    SECTION ("sequence is TTTTGCC")
    {
      DnaString sequence = "TTTTGCC";
      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_01 == actual_bitset);
    }

    SECTION ("sequence is TTTTTGC")
    {
      DnaString sequence = "TTTTTGC";
      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_11 == actual_bitset);
    }

    SECTION ("sequence is TTTTTGCA. The A should not be ignored")
    {
      DnaString sequence = "TTTTTGCA";
      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_10 == actual_bitset);
    }
  }

  SECTION ("Add intron 2 and align using 2 mers")
  {
    int kmer_size = 2;
    String<TVertexDescriptor> order;
    topologicalSort(order, graph);
    TKmerMap kmer_map = kmerifyGraph(order, graph, vertex_vector, free_nodes, edge_ids, kmer_size);

    boost::dynamic_bitset<> expected_bitset_11(num_ids, 3ul);
    boost::dynamic_bitset<> expected_bitset_10(num_ids, 2ul);
    boost::dynamic_bitset<> expected_bitset_01(num_ids, 1ul);
    boost::dynamic_bitset<> expected_bitset_00(num_ids, 0ul);

    REQUIRE(kmer_map.size() == 5);

    DnaString kmer1 = "TT";
    REQUIRE(kmer_map.count(kmer1) == 1);
    REQUIRE(kmer_map[kmer1].size() == 4);
    REQUIRE(kmer_map[kmer1][0].start_vertex == 8);
    REQUIRE(kmer_map[kmer1][0].end_vertex == 9);
    REQUIRE(kmer_map[kmer1][0].id_bits == expected_bitset_11);
    REQUIRE(kmer_map[kmer1][1].start_vertex == 9);
    REQUIRE(kmer_map[kmer1][1].end_vertex == 10);
    REQUIRE(kmer_map[kmer1][1].id_bits == expected_bitset_11);
    REQUIRE(kmer_map[kmer1][2].start_vertex == 10);
    REQUIRE(kmer_map[kmer1][2].end_vertex == 11);
    REQUIRE(kmer_map[kmer1][2].id_bits == expected_bitset_11);
    REQUIRE(kmer_map[kmer1][3].start_vertex == 11);
    REQUIRE(kmer_map[kmer1][3].end_vertex == 2);
    REQUIRE(kmer_map[kmer1][3].id_bits == expected_bitset_11);

    DnaString kmer2 = "TG";
    REQUIRE(kmer_map.count(kmer2) == 1);
    REQUIRE(kmer_map[kmer2].size() == 1);
    REQUIRE(kmer_map[kmer2][0].start_vertex == 2);
    REQUIRE(kmer_map[kmer2][0].end_vertex == 3);
    REQUIRE(kmer_map[kmer2][0].id_bits == expected_bitset_11);

    DnaString kmer3 = "GC";
    REQUIRE(kmer_map.count(kmer3) == 1);
    REQUIRE(kmer_map[kmer3].size() == 2);
    REQUIRE(kmer_map[kmer3][0].start_vertex == 3);
    REQUIRE(kmer_map[kmer3][0].end_vertex == 5);
    REQUIRE(kmer_map[kmer3][0].id_bits == expected_bitset_10);
    REQUIRE(kmer_map[kmer3][1].start_vertex == 3);
    REQUIRE(kmer_map[kmer3][1].end_vertex == 4);
    REQUIRE(kmer_map[kmer3][1].id_bits == expected_bitset_01);

    DnaString kmer4 = "CC";
    REQUIRE(kmer_map.count(kmer4) == 1);
    REQUIRE(kmer_map[kmer4].size() == 1);
    REQUIRE(kmer_map[kmer4][0].start_vertex == 4);
    REQUIRE(kmer_map[kmer4][0].end_vertex == 5);
    REQUIRE(kmer_map[kmer4][0].id_bits == expected_bitset_01);

    DnaString kmer5 = "CA";
    REQUIRE(kmer_map.count(kmer5) == 1);
    REQUIRE(kmer_map[kmer5].size() == 1);
    REQUIRE(kmer_map[kmer5][0].start_vertex == 5);
    REQUIRE(kmer_map[kmer5][0].end_vertex == 6);
    REQUIRE(kmer_map[kmer5][0].id_bits == expected_bitset_10);


    SECTION ("sequence is TTTTTGCC")
    {
      DnaString sequence = "TTTTTGCC";

      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_01 == actual_bitset);
    }

    SECTION ("sequence is TGCC")
    {
      DnaString sequence = "TGCC";

      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_01 == actual_bitset);
    }

    SECTION ("sequence is TTGC")
    {
      DnaString sequence = "TTGC";

      std::vector<TVertexDescriptor> matching_vertices;
      boost::dynamic_bitset<> actual_bitset =
        align_kmer_to_graph (sequence,
                             num_ids,
                             kmer_map,
                             vertex_vector,
                             0,
                             kmer_size
                            );

      REQUIRE (expected_bitset_11 == actual_bitset);
    }
  }
}


TEST_CASE ("Try to align with backward alignments")
{
  std::ostringstream ss;
  ss << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case6.fa";
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

  std::size_t num_ids = ids.size();

  boost::dynamic_bitset<> expected_bitset_11(num_ids, 3ul);
  boost::dynamic_bitset<> expected_bitset_10(num_ids, 2ul);
  boost::dynamic_bitset<> expected_bitset_01(num_ids, 1ul);
  boost::dynamic_bitset<> expected_bitset_00(num_ids, 0ul);

  int kmer_size = 3;
  String<TVertexDescriptor> order;
  topologicalSort(order, graph);

  TKmerMap kmer_map = kmerifyGraph(order, graph, vertex_vector, free_nodes, edge_ids, kmer_size);


  SECTION ("sequence is TTG")
  {
    DnaString sequence = "TTG";
    std::vector<TVertexDescriptor> matching_vertices;
    boost::dynamic_bitset<> actual_bitset =
      align_kmer_to_graph (sequence,
                           num_ids,
                           kmer_map,
                           vertex_vector,
                           0,
                           kmer_size
                          );

    REQUIRE (expected_bitset_11 == actual_bitset);
  }

  SECTION ("sequence is TTG")
  {
    DnaString sequence = "TTG";
    std::vector<TVertexDescriptor> matching_vertices;
    boost::dynamic_bitset<> actual_bitset =
      align_kmer_to_graph (sequence,
                           num_ids,
                           kmer_map,
                           vertex_vector,
                           0,
                           kmer_size
                          );

    REQUIRE (expected_bitset_11 == actual_bitset);
  }

  SECTION ("sequence is TTGGG")
  {
    DnaString sequence = "TTGGG";
    std::vector<TVertexDescriptor> matching_vertices;
    boost::dynamic_bitset<> actual_bitset =
      align_kmer_to_graph (sequence,
                           num_ids,
                           kmer_map,
                           vertex_vector,
                           2,
                           kmer_size
                          );

    REQUIRE (expected_bitset_11 == actual_bitset);
  }

  SECTION ("sequence is TTTGGG, start at index 2")
  {
    DnaString sequence = "TTTGGG";
    std::vector<TVertexDescriptor> matching_vertices;
    boost::dynamic_bitset<> actual_bitset =
      align_kmer_to_graph (sequence,
                           num_ids,
                           kmer_map,
                           vertex_vector,
                           2,
                           kmer_size
                          );

    REQUIRE (expected_bitset_10 == actual_bitset);
  }

  SECTION ("sequence is TTTGGG, start at index 1")
  {
    DnaString sequence = "TTTGGG";
    std::vector<TVertexDescriptor> matching_vertices;
    boost::dynamic_bitset<> actual_bitset =
      align_kmer_to_graph (sequence,
                           num_ids,
                           kmer_map,
                           vertex_vector,
                           1,
                           kmer_size
                          );

    REQUIRE (expected_bitset_10 == actual_bitset);
  }
}



TEST_CASE("find_best_kmer should find the kmer with the best quality")
{
  // !"#$%&'()*+.-./0123456789:;<=>?@ABCDEFGHIJ
  unsigned k = 4;

  SECTION ("Last kmer has the highest quality")
  {
    String<char> qual = "44440";
    REQUIRE(0 == find_best_kmer(qual, k));
  }
  
  SECTION ("First kmer has the highest quality")
  {
    String<char> qual = "04444";
    REQUIRE(1 == find_best_kmer(qual, k));
  }

  SECTION ("Middle kmer has the highest quality")
  {
    String<char> qual = "044440";
    REQUIRE(1 == find_best_kmer(qual, k));
  }

  SECTION ("Just another test case")
  {
    String<char> qual = "0000099990000";
    REQUIRE(5 == find_best_kmer(qual, k));
  }

  SECTION ("Don't break if with highest quality with k = 16")
  {
    k = 16;
    String<char> qual = "0000JJJJJJJJJJJJJJJJ0000";
    REQUIRE(4 == find_best_kmer(qual, k));
  }

  SECTION ("*,A,A,A,A is far worse than 1,1,1,1,1 even though its sum is higher")
  {
    k = 5;
    String<char> qual = "000033333000000*AAAA";
    REQUIRE(4 == find_best_kmer(qual, k));
  }
}

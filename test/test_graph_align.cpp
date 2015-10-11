#include <stdio.h>
#include <climits>

#include "../src/graph_builder.hpp"
#include "../src/graph_align.hpp"
#include "../src/graph_io.hpp"
#include "constants.hpp"
#include "catch.hpp"


using namespace seqan;


/*
 * initializeScoreMatrixAndBacktracker
 **/
TEST_CASE ("initializeExactScoreMatrixAndBacktracker should initialize the score matrix and the backtracker with a simple case")
{
  std::size_t seq_length = 1;
  std::size_t number_of_nodes = 1;
  std::vector<ExactBacktracker> backtracker;

  SECTION ("backtracker requirements")
  {
    initializeExactScoreMatrixAndBacktracker(seq_length, number_of_nodes, backtracker);

    REQUIRE (backtracker.size() == 1);
    REQUIRE (backtracker[0].nodes.size() == 1);
    REQUIRE (backtracker[0].nodes[0] == 0);
    REQUIRE (backtracker[0].match.size() == 1);
    REQUIRE (backtracker[0].match[0] == false);
  }

  seq_length = 2;
  number_of_nodes = 4;

  SECTION ("exact backtracker with a more complicated case")
  {
    initializeExactScoreMatrixAndBacktracker(seq_length, number_of_nodes, backtracker);

    REQUIRE (backtracker.size() == 4);
    REQUIRE (backtracker[0].nodes.size() == 2);
    REQUIRE (backtracker[0].match.size() == 2);

    for (TVertexDescriptor node = 0 ; node < 4 ; ++node)
    {
      REQUIRE (backtracker[node].nodes[0] == node);
      REQUIRE (backtracker[node].nodes[1] == node);
      REQUIRE (backtracker[node].match[0] == false);
      REQUIRE (backtracker[node].match[1] == false);
    }
  }
}


/*
 * getScoreVector with exact matches
 **/
TEST_CASE("getScoreVector with exact matches and a single reference")
{
  DnaString sequence = "GATA";
  std::size_t seq_length = length(sequence);
  std::size_t number_of_nodes = 3;
  std::vector< ExactBacktracker > backtracker;

  initializeExactScoreMatrixAndBacktracker(seq_length, number_of_nodes, backtracker);

  SECTION ("Reference as G")
  {
    TVertexDescriptor source_vertex = 0;
    TVertexDescriptor target_vertex = 1;

    Dna reference = Dna('G');

    bool done =
    getScoreVector(backtracker[source_vertex],
                   backtracker[target_vertex],
                   sequence,
                   reference,
                   source_vertex
                  );

    REQUIRE (!done);
    REQUIRE (backtracker[target_vertex].nodes[0] == source_vertex);
    REQUIRE (backtracker[target_vertex].nodes[1] == target_vertex);
    REQUIRE (backtracker[target_vertex].nodes[2] == target_vertex);
    REQUIRE (backtracker[target_vertex].nodes[3] == target_vertex);

    REQUIRE (backtracker[target_vertex].match[0] == true);
    REQUIRE (backtracker[target_vertex].match[1] == false);
    REQUIRE (backtracker[target_vertex].match[2] == false);
    REQUIRE (backtracker[target_vertex].match[3] == false);

    SECTION ("G and then A")
    {
      reference = Dna('A');
      source_vertex = 1;
      target_vertex = 2;

      bool done =
      getScoreVector(backtracker[source_vertex],
                     backtracker[target_vertex],
                     sequence,
                     reference,
                     source_vertex
                    );

      REQUIRE (!done);
      REQUIRE (backtracker[target_vertex].nodes[0] == target_vertex);
      REQUIRE (backtracker[target_vertex].nodes[1] == source_vertex);
      REQUIRE (backtracker[target_vertex].nodes[2] == target_vertex);
      REQUIRE (backtracker[target_vertex].nodes[3] == target_vertex);

      REQUIRE (backtracker[target_vertex].match[0] == false);
      REQUIRE (backtracker[target_vertex].match[1] == true);
      REQUIRE (backtracker[target_vertex].match[2] == false);
      REQUIRE (backtracker[target_vertex].match[3] == false);
    }

    SECTION ("G and then T")
    {
      reference = Dna('T');
      source_vertex = 1;
      target_vertex = 2;

      bool done =
      getScoreVector(backtracker[source_vertex],
                     backtracker[target_vertex],
                     sequence,
                     reference,
                     source_vertex
                    );

      REQUIRE (!done);
      REQUIRE (backtracker[target_vertex].nodes[0] == target_vertex);
      REQUIRE (backtracker[target_vertex].nodes[1] == target_vertex);
      REQUIRE (backtracker[target_vertex].nodes[2] == target_vertex);
      REQUIRE (backtracker[target_vertex].nodes[3] == target_vertex);

      REQUIRE (backtracker[target_vertex].match[0] == false);
      REQUIRE (backtracker[target_vertex].match[1] == false);
      REQUIRE (backtracker[target_vertex].match[2] == false);
      REQUIRE (backtracker[target_vertex].match[3] == false);
    }
  }

  SECTION ("Reference as T")
  {
    TVertexDescriptor source_vertex = 0;
    TVertexDescriptor target_vertex = 1;

    Dna reference = Dna('T');

    getScoreVector(backtracker[source_vertex],
                   backtracker[target_vertex],
                   sequence,
                   reference,
                   source_vertex
                  );

    REQUIRE (backtracker[target_vertex].nodes[0] == target_vertex);
    REQUIRE (backtracker[target_vertex].nodes[1] == target_vertex);
    REQUIRE (backtracker[target_vertex].nodes[2] == target_vertex);
    REQUIRE (backtracker[target_vertex].nodes[3] == target_vertex);

    REQUIRE (backtracker[target_vertex].match[0] == false);
    REQUIRE (backtracker[target_vertex].match[1] == false);
    REQUIRE (backtracker[target_vertex].match[2] == false);
    REQUIRE (backtracker[target_vertex].match[3] == false);
  }

}


/*
 * alignToGraphExact
 **/
TEST_CASE ("alignToGraphExact should return a backtracker ")
{
  SECTION ("Test case 1")
  {
    // Graph creation
    std::ostringstream ss;
    ss << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case1.fa";
    std::string tmp_str = ss.str();
    const char* alignment_file = tmp_str.c_str();

    std::vector<std::string> ids;
    std::vector<VertexLabels> vertex_vector;
    boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

    TVertexDescriptor begin_vertex;
    TGraph graph = createGraph(alignment_file, vertex_vector, edge_ids, ids, begin_vertex);

    String<TVertexDescriptor> order;
    topologicalSort(order, graph);
    std::vector< ExactBacktracker > backtracker;


    DnaString sequence = "GATA";

    initializeExactScoreMatrixAndBacktracker(length(sequence), length(order), backtracker);

    boost::unordered_set<TVertexDescriptor> free_nodes;
    free_nodes.insert(1);

    std::vector<TVertexDescriptor> matching_vertices;
    alignToGraphExact (sequence,
                       order,
                       graph,
                       matching_vertices,
                       vertex_vector,
                       backtracker,
                       free_nodes
                      );

    REQUIRE(length(order) == 4);
    REQUIRE(length(sequence) == 4);
    REQUIRE(matching_vertices.size() == 0);

    REQUIRE(backtracker.size() == 4);
    REQUIRE(backtracker[0].match[0] == false);
    REQUIRE(backtracker[0].match[1] == false);
    REQUIRE(backtracker[0].match[2] == false);
    REQUIRE(backtracker[0].match[3] == false);

    REQUIRE(backtracker[1].match[0] == false);
    REQUIRE(backtracker[1].match[1] == false);
    REQUIRE(backtracker[1].match[2] == false);
    REQUIRE(backtracker[1].match[3] == false);

    REQUIRE(backtracker[2].match[0] == false);
    REQUIRE(backtracker[2].match[1] == false);
    REQUIRE(backtracker[2].match[2] == false);
    REQUIRE(backtracker[2].match[3] == false);

    REQUIRE(backtracker[3].match[0] == false);
    REQUIRE(backtracker[3].match[1] == false);
    REQUIRE(backtracker[3].match[2] == false);
    REQUIRE(backtracker[3].match[3] == false);
  }

  SECTION ("Test case 2")
  {
    DnaString sequence = "GATA";

    std::ostringstream ss;
    ss << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case2.fa";
    std::string tmp_str = ss.str();
    const char* alignment_file = tmp_str.c_str();

    std::vector<std::string> ids;
    std::vector<VertexLabels> vertex_vector;
    boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

    TVertexDescriptor begin_vertex;
    TGraph graph = createGraph(alignment_file, vertex_vector, edge_ids, ids, begin_vertex);

    String<TVertexDescriptor> order;
    topologicalSort(order, graph);
    std::vector< ExactBacktracker > backtracker;

    initializeExactScoreMatrixAndBacktracker(length(sequence), length(order), backtracker);

    boost::unordered_set<TVertexDescriptor> free_nodes;
    free_nodes.insert(1);

    std::vector<TVertexDescriptor> matching_vertices;
    alignToGraphExact (sequence,
                       order,
                       graph,
                       matching_vertices,
                       vertex_vector,
                       backtracker,
                       free_nodes
                      );

    REQUIRE(length(order) == 6);
    REQUIRE(length(sequence) == 4);

    REQUIRE(backtracker.size() == 6);
    REQUIRE(backtracker[0].match[0] == false);
    REQUIRE(backtracker[0].match[1] == false);
    REQUIRE(backtracker[0].match[2] == false);
    REQUIRE(backtracker[0].match[3] == false);

    REQUIRE(backtracker[1].match[0] == false);
    REQUIRE(backtracker[1].match[1] == false);
    REQUIRE(backtracker[1].match[2] == false);
    REQUIRE(backtracker[1].match[3] == true);

    REQUIRE(backtracker[2].match[0] == true);
    REQUIRE(backtracker[2].match[1] == false);
    REQUIRE(backtracker[2].match[2] == false);
    REQUIRE(backtracker[2].match[3] == false);

    REQUIRE(backtracker[3].match[0] == false);
    REQUIRE(backtracker[3].match[1] == true);
    REQUIRE(backtracker[3].match[2] == false);
    REQUIRE(backtracker[3].match[3] == false);

    REQUIRE(backtracker[4].match[0] == false);
    REQUIRE(backtracker[4].match[1] == false);
    REQUIRE(backtracker[4].match[2] == true);
    REQUIRE(backtracker[4].match[3] == false);

    REQUIRE(backtracker[5].match[0] == false);
    REQUIRE(backtracker[5].match[1] == false);
    REQUIRE(backtracker[5].match[2] == false);
    REQUIRE(backtracker[5].match[3] == true);

    REQUIRE(matching_vertices.size() == 1);
    REQUIRE(matching_vertices[0] == 5);
  }

  SECTION ("Test case 3")
  {
    DnaString sequence = "GATA";

    std::ostringstream ss;
    ss << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case3.fa";
    std::string tmp_str = ss.str();
    const char* alignment_file = tmp_str.c_str();

    std::vector<std::string> ids;
    std::vector<VertexLabels> vertex_vector;
    boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

    TVertexDescriptor begin_vertex;
    TGraph graph = createGraph(alignment_file, vertex_vector, edge_ids, ids, begin_vertex);

    String<TVertexDescriptor> order;
    topologicalSort(order, graph);
    std::vector< ExactBacktracker > backtracker;

    initializeExactScoreMatrixAndBacktracker(length(sequence), length(order), backtracker);

    boost::unordered_set<TVertexDescriptor> free_nodes;
    free_nodes.insert(1);

    std::vector<TVertexDescriptor> matching_vertices;
    alignToGraphExact (sequence,
                       order,
                       graph,
                       matching_vertices,
                       vertex_vector,
                       backtracker,
                       free_nodes
                      );

    REQUIRE(length(order) == 6);
    REQUIRE(length(sequence) == 4);

    REQUIRE(backtracker.size() == 6);
    REQUIRE(backtracker[0].match[0] == false);
    REQUIRE(backtracker[0].match[1] == false);
    REQUIRE(backtracker[0].match[2] == false);
    REQUIRE(backtracker[0].match[3] == false);

    REQUIRE(backtracker[1].match[0] == false);
    REQUIRE(backtracker[1].match[1] == false);
    REQUIRE(backtracker[1].match[2] == true);
    REQUIRE(backtracker[1].match[3] == true);

    REQUIRE(backtracker[2].match[0] == true);
    REQUIRE(backtracker[2].match[1] == false);
    REQUIRE(backtracker[2].match[2] == false);
    REQUIRE(backtracker[2].match[3] == false);

    REQUIRE(backtracker[3].match[0] == false);
    REQUIRE(backtracker[3].match[1] == true);
    REQUIRE(backtracker[3].match[2] == false);
    REQUIRE(backtracker[3].match[3] == false);

    REQUIRE(backtracker[4].match[0] == false);
    REQUIRE(backtracker[4].match[1] == false);
    REQUIRE(backtracker[4].match[2] == true);
    REQUIRE(backtracker[4].match[3] == false);

    REQUIRE(backtracker[5].match[0] == false);
    REQUIRE(backtracker[5].match[1] == false);
    REQUIRE(backtracker[5].match[2] == false);
    REQUIRE(backtracker[5].match[3] == true);

    REQUIRE(matching_vertices.size() == 1);
    REQUIRE(matching_vertices[0] == 5);
  }
}


/*
 * alignToGraph
 **/
TEST_CASE ("test case 3 wih exact matches and alignGraph and extendGraph")
{
  std::ostringstream ss;
  ss << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case4_e3.fa";
  std::string tmp_str = ss.str();
  const char* alignment_file_e3 = tmp_str.c_str();

  std::vector<std::string> ids;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;
  boost::unordered_set<TVertexDescriptor> free_nodes;
  free_nodes.insert(1);

  TVertexDescriptor begin_vertex;
  TGraph graph = createGraph(alignment_file_e3, vertex_vector, edge_ids, ids, begin_vertex);

  std::ostringstream ss2;
  ss2 << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case4_i2.fa";
  std::string tmp_str2 = ss2.str();
  const char* alignment_file_i2 = tmp_str2.c_str();
  TVertexDescriptor new_begin_vertex;
  free_nodes.insert(begin_vertex);
  extendGraph(graph, alignment_file_i2, vertex_vector, new_begin_vertex, begin_vertex);
  std::vector<ExactBacktracker> backtracker;

  SECTION ("exon 3 and intron 2")
  {
    String<TVertexDescriptor> order;
    topologicalSort(order, graph);

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
  }

  
  std::ostringstream ss3;
  ss3 << gyper_SOURCE_DIRECTORY << "/test/alignments/test_case4_e2.fa";
  std::string tmp_str3 = ss3.str();
  const char* alignment_file_e2 = tmp_str3.c_str();
  free_nodes.insert(begin_vertex);
  extendGraph(graph, alignment_file_e2, vertex_vector, edge_ids, new_begin_vertex, begin_vertex);

  SECTION ("exon 2 and 3, and intron 2")
  {
    String<TVertexDescriptor> order;
    topologicalSort(order, graph);
    REQUIRE (length(order) == 16);

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

      REQUIRE (backtracker.size() == length(order));
      REQUIRE (backtracker[0].match.size() == length(sequence));

      REQUIRE (backtracker[0].match[0] == false);
      REQUIRE (backtracker[0].match[1] == false);
      REQUIRE (backtracker[0].match[2] == false);
      REQUIRE (backtracker[0].match[3] == false);

      REQUIRE (backtracker[1].match[0] == false);
      REQUIRE (backtracker[1].match[1] == false);
      REQUIRE (backtracker[1].match[2] == true);
      REQUIRE (backtracker[1].match[3] == true);

      REQUIRE (backtracker[2].match[0] == true);
      REQUIRE (backtracker[2].match[1] == false);
      REQUIRE (backtracker[2].match[2] == false);
      REQUIRE (backtracker[2].match[3] == false);

      REQUIRE (backtracker[3].match[0] == false);
      REQUIRE (backtracker[3].match[1] == true);
      REQUIRE (backtracker[3].match[2] == false);
      REQUIRE (backtracker[3].match[3] == false);

      REQUIRE (backtracker[4].match[0] == false);
      REQUIRE (backtracker[4].match[1] == false);
      REQUIRE (backtracker[4].match[2] == true);
      REQUIRE (backtracker[4].match[3] == false);

      REQUIRE (backtracker[5].match[0] == false);
      REQUIRE (backtracker[5].match[1] == false);
      REQUIRE (backtracker[5].match[2] == false);
      REQUIRE (backtracker[5].match[3] == true);

      boost::dynamic_bitset<> perfect_match =
        backTrackAndCount(ids, backtracker, matching_vertices[0], edge_ids);

      boost::dynamic_bitset<> expected(2, 1ul);
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

      REQUIRE (backtracker[0].match[0] == true);
      REQUIRE (backtracker[0].match[1] == false);
      REQUIRE (backtracker[0].match[2] == false);
      REQUIRE (backtracker[0].match[3] == false);
      REQUIRE (backtracker[0].match[4] == false);

      REQUIRE (backtracker[1].match[0] == true);
      REQUIRE (backtracker[1].match[1] == false);
      REQUIRE (backtracker[1].match[2] == false);
      REQUIRE (backtracker[1].match[3] == true);

      REQUIRE (backtracker[2].match[1] == true);
      REQUIRE (backtracker[3].match[2] == true);
      REQUIRE (backtracker[4].match[3] == true);
      REQUIRE (backtracker[5].match[4] == true);

      boost::dynamic_bitset<> perfect_match =
        backTrackAndCount(ids, backtracker, matching_vertices[0], edge_ids);

      boost::dynamic_bitset<> expected(2, 1ul);
      REQUIRE (expected == perfect_match);
    }

    SECTION ("sequence is CTTAGATA")
    {
      DnaString sequence = "CTTAGATA";

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

      REQUIRE (backtracker[14].match[0] == true);
      REQUIRE (backtracker[6].match[0] == true);
      REQUIRE (backtracker[6].match[1] == false);
      REQUIRE (backtracker[6].match[2] == false);

      REQUIRE (backtracker[0].match[0] == false);
      REQUIRE (backtracker[0].match[1] == false);
      REQUIRE (backtracker[0].match[2] == false);
      REQUIRE (backtracker[0].match[3] == true);
      REQUIRE (backtracker[0].match[4] == false);
      REQUIRE (backtracker[0].match[5] == false);
      REQUIRE (backtracker[0].match[6] == false);
      REQUIRE (backtracker[0].match[7] == false);

      REQUIRE (backtracker[2].match[4] == true);
      REQUIRE (backtracker[3].match[5] == true);
      REQUIRE (backtracker[4].match[6] == true);
      REQUIRE (backtracker[5].match[7] == true);

      boost::dynamic_bitset<> perfect_match =
        backTrackAndCount(ids, backtracker, matching_vertices[0], edge_ids);

      boost::dynamic_bitset<> expected(2, 1ul);
      REQUIRE (expected == perfect_match);
    }

    SECTION ("If sequence is ATTAGATA, there shouldn't be any matches")
    {
      DnaString sequence = "ATTAGATA";

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

      REQUIRE (matching_vertices.size() == 0);

      REQUIRE (backtracker[14].match[0] == false);
      REQUIRE (backtracker[6].match[0] == false);
      REQUIRE (backtracker[6].match[1] == false);
      REQUIRE (backtracker[6].match[2] == false);

      REQUIRE (backtracker[9].match[0] == true);
      REQUIRE (backtracker[0].match[0] == true);
      REQUIRE (backtracker[0].match[1] == false);
      REQUIRE (backtracker[0].match[2] == false);
      REQUIRE (backtracker[0].match[3] == false);
      REQUIRE (backtracker[0].match[4] == false);
      REQUIRE (backtracker[0].match[5] == false);
      REQUIRE (backtracker[0].match[6] == false);
      REQUIRE (backtracker[0].match[7] == false);

      REQUIRE (backtracker[2].match[4] == false);
      REQUIRE (backtracker[3].match[5] == false);
      REQUIRE (backtracker[4].match[6] == false);
      REQUIRE (backtracker[5].match[7] == false);

      REQUIRE (backtracker[3].match[0] == true);
      REQUIRE (backtracker[4].match[1] == true);
      REQUIRE (backtracker[1].match[1] == true);

      if (matching_vertices.size() != 0)
      {
        backTrackAndCount(ids, backtracker, matching_vertices[0], edge_ids);
      }
    }

  }
}



TEST_CASE ("Vcf file out")
{
  SECTION ("Write VCF file")
  {
    std::string pn = "ABKDWWG";
    std::vector<std::string> ids;
    ids.push_back("HLA-A*01:01");
    ids.push_back("HLA-A*01:06");
    ids.push_back("HLA-A*18:01");
    ids.push_back("HLA-A*102:01");
    std::vector<std::vector<double> > seq_scores;
    seq_scores.resize(ids.size());
    for (unsigned pos = 0 ; pos < ids.size() ; ++pos)
    {
      seq_scores[pos].resize(ids.size());
    }
    // std::cout << seq_scores.size() << " " << seq_scores[0].size() << std::endl;
    
    seq_scores[0][0] = 8.5;
    seq_scores[1][0] = 9.75;
    seq_scores[1][1] = 2;
    seq_scores[2][0] = 10.0;
    seq_scores[2][1] = 5;
    seq_scores[2][2] = 0.0;
    seq_scores[3][0] = 1.5;
    seq_scores[3][1] = 3.0;
    seq_scores[3][2] = 4.5;
    seq_scores[3][3] = 6.0;

    double max_score = 10;
    // std::ostringstream my_ss;
    CharString file_name = "example_test.vcf";
    testVcfFile(file_name, pn, ids, seq_scores, max_score);
  }
}

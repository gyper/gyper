#include <stdio.h>
#include <climits>

#include "../src/graph_builder.hpp"
#include "../src/graph_align.hpp"
#include "constants.hpp"
#include "catch.hpp"


using namespace seqan;


TEST_CASE ("addSequenceToGraph should be able to add a single sequence")
{
  TGraph g;
  
  SECTION ("graph should be empty to begin with")
  {
    REQUIRE (numVertices(g) == 0);
    REQUIRE (numEdges(g) == 0);
  }

  CharString sequence = "TAG";
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  std::vector<VertexLabels> vertex_vector;
  TVertexDescriptor begin_vertex = 0;

  addInitialAndEndVertex(g, vertex_label_map, vertex_vector, begin_vertex);
  addSequenceToGraph(g, 0, 1, sequence, vertex_label_map, vertex_vector, edge_ids);

  SECTION ("graph with sequence 'TAG' should now have 4 vertices and 3 edges")
  {
    REQUIRE (numVertices(g) == 5);
    REQUIRE (numEdges(g) == 4);
  }
}

void
changePair(std::pair<TVertexDescriptor, TVertexDescriptor> & my_pair, TVertexDescriptor new_first, TVertexDescriptor new_second)
{
  my_pair.first = new_first;
  my_pair.second = new_second;
}

TEST_CASE ("addSequenceToGraph should be able to add a single sequence with gaps")
{
  TGraph g;
  CharString sequence = "TA-G";
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  std::vector<String<CharString> > id_vector;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  TVertexDescriptor begin_vertex = 0;
  addInitialAndEndVertex(g, vertex_label_map, vertex_vector, begin_vertex);
  addSequenceToGraph(g, 0, 1, sequence, vertex_label_map, vertex_vector, edge_ids);

  String<bool> matrix;
  getAdjacencyMatrix(g, matrix);

  SECTION ("graph with sequence 'TA-G' should now have 4 vertices and 3 edges")
  {
    REQUIRE (numVertices(g) == 5);
    REQUIRE (numEdges(g) == 4);
  }

  SECTION ("there should be an edge between the labels")
  {
    REQUIRE (!hasEdge(matrix, 0, 1, 5));
    REQUIRE (hasEdge(matrix, 0, 2, 5));
    REQUIRE (hasEdge(matrix, 2, 3, 5));
    REQUIRE (hasEdge(matrix, 3, 4, 5));
    REQUIRE (hasEdge(matrix, 4, 1, 5));

    REQUIRE (!hasEdge(matrix, 0, 1));
    REQUIRE (hasEdge(matrix, 0, 2));
    REQUIRE (hasEdge(matrix, 2, 3));
    REQUIRE (hasEdge(matrix, 3, 4));
    REQUIRE (hasEdge(matrix, 4, 1));
  }

  SECTION ("there shouldn't be any more edges in that matrix")
  {
    unsigned sum = 0;
    for (unsigned i = 0 ; i < length(matrix) ; ++i)
    {
      if ( matrix[i] != 0 )
        ++sum;
    }
    REQUIRE (sum == numEdges(g));
    REQUIRE (sum == 4);
  }

  SECTION ("edge_ids should have the added id")
  {
    REQUIRE(edge_ids.size() == 4);
    std::pair<TVertexDescriptor, TVertexDescriptor> my_pair(0, 2);
    REQUIRE (edge_ids.count(my_pair) == 1);
    REQUIRE (edge_ids[my_pair].to_ulong() == 1);

    changePair(my_pair, 2, 3);
    REQUIRE (edge_ids.count(my_pair) == 1);
    REQUIRE (edge_ids[my_pair].to_ulong() == 1);

    changePair(my_pair, 3, 4);
    REQUIRE (edge_ids.count(my_pair) == 1);
    REQUIRE (edge_ids[my_pair].to_ulong() == 1);

    changePair(my_pair, 4, 1);
    REQUIRE (edge_ids.count(my_pair) == 1);
    REQUIRE (edge_ids[my_pair].to_ulong() == 1);
  }
}


TEST_CASE ("addSequenceToGraph should be able to add more than one sequence")
{
  TGraph g;
  CharString sequence1 = "TA-G";
  CharString sequence2 = "TAGG";
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  std::vector<String<CharString> > id_vector;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  TVertexDescriptor begin_vertex = 0;
  addInitialAndEndVertex(g, vertex_label_map, vertex_vector, begin_vertex);
  addSequenceToGraph(g, 0, 2, sequence1, vertex_label_map, vertex_vector, edge_ids);
  addSequenceToGraph(g, 1, 2, sequence2, vertex_label_map, vertex_vector, edge_ids);

  SECTION ("graph should now have 2 more edges and 1 more vertix")
  {
    REQUIRE (numVertices(g) == 6);
    REQUIRE (numEdges(g) == 6);
  }

  SECTION ("if we add an existing sequence in again, no new vertices or edges will be added")
  {
    addSequenceToGraph(g, 1, 2, sequence2, vertex_label_map, vertex_vector, edge_ids);
    REQUIRE (numVertices(g) == 6);
    REQUIRE (numEdges(g) == 6);
  }

  String<bool> matrix;
  getAdjacencyMatrix(g, matrix);

  SECTION ("there should be an edge between the labels")
  {
    REQUIRE (length(matrix) == 6*6);

    // Test1: 0 -T> 2 -A> 3 -G> 4 -> 1
    // Test2: 0 -T> 2 -A> 3 -G> 5 -G> 4 -> 1
    REQUIRE (hasEdge(matrix, 0, 0) == false);
    REQUIRE (hasEdge(matrix, 1, 0) == false);
    REQUIRE (hasEdge(matrix, 2, 0) == false);
    REQUIRE (hasEdge(matrix, 3, 0) == false);
    REQUIRE (hasEdge(matrix, 4, 0) == false);
    REQUIRE (hasEdge(matrix, 5, 0) == false);

    REQUIRE (hasEdge(matrix, 0, 1) == false);
    REQUIRE (hasEdge(matrix, 1, 1) == false);
    REQUIRE (hasEdge(matrix, 2, 1) == false);
    REQUIRE (hasEdge(matrix, 3, 1) == false);
    REQUIRE (hasEdge(matrix, 4, 1) == true);
    REQUIRE (hasEdge(matrix, 5, 1) == false);

    REQUIRE (hasEdge(matrix, 0, 2) == true);
    REQUIRE (hasEdge(matrix, 1, 2) == false);
    REQUIRE (hasEdge(matrix, 2, 2) == false);
    REQUIRE (hasEdge(matrix, 3, 2) == false);
    REQUIRE (hasEdge(matrix, 4, 2) == false);
    REQUIRE (hasEdge(matrix, 5, 2) == false);

    REQUIRE (hasEdge(matrix, 0, 3) == false);
    REQUIRE (hasEdge(matrix, 1, 3) == false);
    REQUIRE (hasEdge(matrix, 2, 3) == true);
    REQUIRE (hasEdge(matrix, 3, 3) == false);
    REQUIRE (hasEdge(matrix, 4, 3) == false);
    REQUIRE (hasEdge(matrix, 5, 3) == false);

    REQUIRE (hasEdge(matrix, 0, 4) == false);
    REQUIRE (hasEdge(matrix, 1, 4) == false);
    REQUIRE (hasEdge(matrix, 2, 4) == false);
    REQUIRE (hasEdge(matrix, 3, 4) == true);
    REQUIRE (hasEdge(matrix, 4, 4) == false);
    REQUIRE (hasEdge(matrix, 5, 4) == true);

    REQUIRE (hasEdge(matrix, 0, 5) == false);
    REQUIRE (hasEdge(matrix, 1, 5) == false);
    REQUIRE (hasEdge(matrix, 2, 5) == false);
    REQUIRE (hasEdge(matrix, 3, 5) == true);
    REQUIRE (hasEdge(matrix, 4, 5) == false);
    REQUIRE (hasEdge(matrix, 5, 5) == false);
  }

  SECTION ("there shouldn't be any more edges in that matrix")
  {
    unsigned sum = 0;
    for (unsigned i = 0 ; i < length(matrix) ; ++i)
    {
      if ( matrix[i] != 0 )
        ++sum;
    }
    REQUIRE (sum == 6);
  }
}


TEST_CASE ("addSequenceToGraph should return a map that can point to each of the vertices")
{
  TGraph g;
  CharString sequence = "TA-G";
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  TVertexDescriptor begin_vertex = 0;
  addInitialAndEndVertex(g, vertex_label_map, vertex_vector, begin_vertex);
  addSequenceToGraph(g, 0, 1, sequence, vertex_label_map, vertex_vector, edge_ids);

  VertexLabels label1 = { 0, Dna('A') };
  VertexLabels label2 = { 1, Dna('T') };
  VertexLabels label3 = { 2, Dna('A') };
  VertexLabels label4 = { 4, Dna('G') };


  SECTION ("the map should have 4 entries")
  {
    REQUIRE (vertex_label_map.size() == 5);
  }

  SECTION ("there should be a entry that points to the first vertex")
  {
    REQUIRE (vertex_label_map.count(label1) == 1);    
  }

  SECTION ("there should be a entry that points to the second vertex")
  { 
    REQUIRE (vertex_label_map.count(label2) == 1);
  }

  SECTION ("there should be a entry that points to the third vertex")
  { 
    REQUIRE (vertex_label_map.count(label3) == 1);
  }

  SECTION ("there should be a entry that points to the fourth vertex")
  { 
    REQUIRE (vertex_label_map.count(label4) == 1);
  }
}


TEST_CASE ("createVertexMap should return a vector with vertices as keys and vertex labels as value")
{
  TGraph g;
  CharString sequence1 = "TA-G";
  CharString sequence2 = "TAGG";
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  TVertexDescriptor begin_vertex = 0;
  addInitialAndEndVertex(g, vertex_label_map, vertex_vector, begin_vertex);
  addSequenceToGraph(g, 0, 2, sequence1, vertex_label_map, vertex_vector, edge_ids);
  addSequenceToGraph(g, 1, 2, sequence2, vertex_label_map, vertex_vector, edge_ids);

  SECTION ("The size should be the same as the vertex label map")
  {
    REQUIRE (vertex_vector.size() == vertex_label_map.size());
  }

  SECTION ("The labels should be like the following")
  {
    VertexLabels v_label0 = { 0, Dna('A') };
    TVertexDescriptor vertex0 = vertex_label_map[v_label0];
    
    REQUIRE (vertex_vector[vertex0].level == 0);
    REQUIRE (vertex_vector[vertex0].dna == 0);

    VertexLabels v_label1 = { 1, Dna('T') };
    TVertexDescriptor vertex1 = vertex_label_map[v_label1];
    REQUIRE (vertex_vector[vertex1].level == 1);
    REQUIRE (vertex_vector[vertex1].dna == Dna('T'));

    VertexLabels v_label2 = { 2, Dna('A') };
    TVertexDescriptor vertex2 = vertex_label_map[v_label2];
    REQUIRE (vertex_vector[vertex2].level == 2);
    REQUIRE (vertex_vector[vertex2].dna == Dna('A'));

    VertexLabels v_label3 = { 4, Dna('G') };
    TVertexDescriptor vertex3 = vertex_label_map[v_label3];
    REQUIRE (vertex_vector[vertex3].level == 4);
    REQUIRE (vertex_vector[vertex3].dna == Dna('G'));

    VertexLabels v_label4 = { 3, Dna('G') };
    TVertexDescriptor vertex4 = vertex_label_map[v_label4];
    REQUIRE (vertex_vector[vertex4].level == 3);
    REQUIRE (vertex_vector[vertex4].dna == Dna('G'));
  }
}


TEST_CASE ("addSequenceToGraph should work with matrix of boost's dynamic bits")
{
  TGraph g;
  CharString sequence1 = "TA-G";
  unsigned short bit_id1 = 0;
  CharString sequence2 = "TAGG";
  unsigned short bit_id2 = 1;
  unsigned short bit_n = 2;
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  TVertexDescriptor begin_vertex = 0;
  addInitialAndEndVertex(g, vertex_label_map, vertex_vector, begin_vertex);
  addSequenceToGraph(g, bit_id1, bit_n, sequence1, vertex_label_map, vertex_vector, edge_ids);

  SECTION ("graph should have 5 vertices and 4 edges in total after the first sequence")
  {
    REQUIRE (numVertices(g) == 5);
    REQUIRE (numEdges(g) == 4);
    REQUIRE (vertex_label_map.size() == 5);
    REQUIRE (vertex_vector.size() == 5);
  }

  addSequenceToGraph(g, bit_id2, bit_n, sequence2, vertex_label_map, vertex_vector, edge_ids);

  SECTION ("graph should now have 2 more edges and 1 more vertix")
  {
    REQUIRE (numVertices(g) == 6);
    REQUIRE (numEdges(g) == 6);
  }

  SECTION ("if we add an existing sequence in again, no new vertices or edges will be added")
  {
    addSequenceToGraph(g, bit_id2, bit_n, sequence2, vertex_label_map, vertex_vector, edge_ids);
    REQUIRE (numVertices(g) == 6);
    REQUIRE (numEdges(g) == 6);
  }

  SECTION ("there should be an edge between the labels")
  {
    String<bool> matrix;
    getAdjacencyMatrix(g, matrix);

    REQUIRE (length(matrix) == 6*6);

    // Test1: 0 -T> 2 -A> 3 -G> 4 -> 1
    // Test2: 0 -T> 2 -A> 3 -G> 5 -G> 4 -> 1
    REQUIRE (hasEdge(matrix, 0, 0) == false);
    REQUIRE (hasEdge(matrix, 1, 0) == false);
    REQUIRE (hasEdge(matrix, 2, 0) == false);
    REQUIRE (hasEdge(matrix, 3, 0) == false);
    REQUIRE (hasEdge(matrix, 4, 0) == false);
    REQUIRE (hasEdge(matrix, 5, 0) == false);

    REQUIRE (hasEdge(matrix, 0, 1) == false);
    REQUIRE (hasEdge(matrix, 1, 1) == false);
    REQUIRE (hasEdge(matrix, 2, 1) == false);
    REQUIRE (hasEdge(matrix, 3, 1) == false);
    REQUIRE (hasEdge(matrix, 4, 1) == true);
    REQUIRE (hasEdge(matrix, 5, 1) == false);

    REQUIRE (hasEdge(matrix, 0, 2) == true);
    REQUIRE (hasEdge(matrix, 1, 2) == false);
    REQUIRE (hasEdge(matrix, 2, 2) == false);
    REQUIRE (hasEdge(matrix, 3, 2) == false);
    REQUIRE (hasEdge(matrix, 4, 2) == false);
    REQUIRE (hasEdge(matrix, 5, 2) == false);

    REQUIRE (hasEdge(matrix, 0, 3) == false);
    REQUIRE (hasEdge(matrix, 1, 3) == false);
    REQUIRE (hasEdge(matrix, 2, 3) == true);
    REQUIRE (hasEdge(matrix, 3, 3) == false);
    REQUIRE (hasEdge(matrix, 4, 3) == false);
    REQUIRE (hasEdge(matrix, 5, 3) == false);

    REQUIRE (hasEdge(matrix, 0, 4) == false);
    REQUIRE (hasEdge(matrix, 1, 4) == false);
    REQUIRE (hasEdge(matrix, 2, 4) == false);
    REQUIRE (hasEdge(matrix, 3, 4) == true);
    REQUIRE (hasEdge(matrix, 4, 4) == false);
    REQUIRE (hasEdge(matrix, 5, 4) == true);

    REQUIRE (hasEdge(matrix, 0, 5) == false);
    REQUIRE (hasEdge(matrix, 1, 5) == false);
    REQUIRE (hasEdge(matrix, 2, 5) == false);
    REQUIRE (hasEdge(matrix, 3, 5) == true);
    REQUIRE (hasEdge(matrix, 4, 5) == false);
    REQUIRE (hasEdge(matrix, 5, 5) == false);
  }

  SECTION ("there shouldn't be any more edges in that matrix")
  {
    String<bool> matrix;
    getAdjacencyMatrix(g, matrix);

    unsigned sum = 0;
    for (unsigned i = 0 ; i < length(matrix) ; ++i)
    {
      if ( matrix[i] != 0 )
        ++sum;
    }
    REQUIRE (sum == 6);
  }

  SECTION ("there should be entries in the vertex_label_map")
  {
    VertexLabels label1 = { 0, Dna('A') };
    VertexLabels label2 = { 1, Dna('T') };
    VertexLabels label3 = { 2, Dna('A') };
    VertexLabels label4 = { 4, Dna('G') };
    VertexLabels label5 = { 3, Dna('G') };

    REQUIRE (vertex_label_map.size() == 6);
    REQUIRE (vertex_label_map.count(label1) == 1);    
    REQUIRE (vertex_label_map.count(label2) == 1);
    REQUIRE (vertex_label_map.count(label3) == 1);
    REQUIRE (vertex_label_map.count(label4) == 1);
    REQUIRE (vertex_label_map.count(label5) == 1);
  }

  SECTION ("edge_ids should have the added id")
  {
    std::pair<TVertexDescriptor, TVertexDescriptor> pair1(0, 2);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair2(2, 3);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair3(3, 4);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair4(4, 1);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair5(3, 5);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair6(5, 4);

    REQUIRE (edge_ids.size() == 6);
    REQUIRE (edge_ids.count(pair1) == 1);
    REQUIRE (edge_ids.count(pair2) == 1);
    REQUIRE (edge_ids.count(pair3) == 1);
    REQUIRE (edge_ids.count(pair4) == 1);
    REQUIRE (edge_ids.count(pair5) == 1);
    REQUIRE (edge_ids.count(pair6) == 1);

    REQUIRE (edge_ids[pair1].to_ulong() == 3);
    REQUIRE (edge_ids[pair2].to_ulong() == 3);
    REQUIRE (edge_ids[pair3].to_ulong() == 1);
    REQUIRE (edge_ids[pair4].to_ulong() == 3);
    REQUIRE (edge_ids[pair5].to_ulong() == 2);
    REQUIRE (edge_ids[pair6].to_ulong() == 2);
  }
}

/*
TEST_CASE ("createGraph should, well, create a graph")
{
  const char* alignment_file = (char*)"../data/Alignments/test_case5.fa";
  std::vector<std::string> ids;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;
  TVertexDescriptor begin_vertex;
  TGraph g = createGraph(alignment_file, vertex_vector, edge_ids, ids, begin_vertex);

  SECTION ("graph should now have 2 more edges and 1 more vertix")
  {
    REQUIRE (numVertices(g) == 6);
    REQUIRE (numEdges(g) == 6);
  }

  SECTION ("there should be an edge between the labels")
  {
    String<bool> matrix;
    getAdjacencyMatrix(g, matrix);

    REQUIRE (length(matrix) == 6*6);

    // Test1: 0 -T> 2 -A> 3 -G> 4 -> 1
    // Test2: 0 -T> 2 -A> 3 -G> 5 -G> 4 -> 1
    REQUIRE (hasEdge(matrix, 0, 0) == false);
    REQUIRE (hasEdge(matrix, 1, 0) == false);
    REQUIRE (hasEdge(matrix, 2, 0) == false);
    REQUIRE (hasEdge(matrix, 3, 0) == false);
    REQUIRE (hasEdge(matrix, 4, 0) == false);
    REQUIRE (hasEdge(matrix, 5, 0) == false);

    REQUIRE (hasEdge(matrix, 0, 1) == false);
    REQUIRE (hasEdge(matrix, 1, 1) == false);
    REQUIRE (hasEdge(matrix, 2, 1) == false);
    REQUIRE (hasEdge(matrix, 3, 1) == false);
    REQUIRE (hasEdge(matrix, 4, 1) == true);
    REQUIRE (hasEdge(matrix, 5, 1) == false);

    REQUIRE (hasEdge(matrix, 0, 2) == true);
    REQUIRE (hasEdge(matrix, 1, 2) == false);
    REQUIRE (hasEdge(matrix, 2, 2) == false);
    REQUIRE (hasEdge(matrix, 3, 2) == false);
    REQUIRE (hasEdge(matrix, 4, 2) == false);
    REQUIRE (hasEdge(matrix, 5, 2) == false);

    REQUIRE (hasEdge(matrix, 0, 3) == false);
    REQUIRE (hasEdge(matrix, 1, 3) == false);
    REQUIRE (hasEdge(matrix, 2, 3) == true);
    REQUIRE (hasEdge(matrix, 3, 3) == false);
    REQUIRE (hasEdge(matrix, 4, 3) == false);
    REQUIRE (hasEdge(matrix, 5, 3) == false);

    REQUIRE (hasEdge(matrix, 0, 4) == false);
    REQUIRE (hasEdge(matrix, 1, 4) == false);
    REQUIRE (hasEdge(matrix, 2, 4) == false);
    REQUIRE (hasEdge(matrix, 3, 4) == true);
    REQUIRE (hasEdge(matrix, 4, 4) == false);
    REQUIRE (hasEdge(matrix, 5, 4) == true);

    REQUIRE (hasEdge(matrix, 0, 5) == false);
    REQUIRE (hasEdge(matrix, 1, 5) == false);
    REQUIRE (hasEdge(matrix, 2, 5) == false);
    REQUIRE (hasEdge(matrix, 3, 5) == true);
    REQUIRE (hasEdge(matrix, 4, 5) == false);
    REQUIRE (hasEdge(matrix, 5, 5) == false);
  }

  SECTION ("there shouldn't be any more edges in that matrix")
  {
    String<bool> matrix;
    getAdjacencyMatrix(g, matrix);

    unsigned sum = 0;
    for (unsigned i = 0 ; i < length(matrix) ; ++i)
    {
      if ( matrix[i] != 0 )
        ++sum;
    }
    REQUIRE (sum == 6);
  }

  SECTION ("edge_ids should have the added id")
  {
    std::pair<TVertexDescriptor, TVertexDescriptor> pair1(0, 2);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair2(2, 3);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair3(3, 4);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair4(4, 1);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair5(3, 5);
    std::pair<TVertexDescriptor, TVertexDescriptor> pair6(5, 4);

    REQUIRE (edge_ids.size() == 6);
    REQUIRE (edge_ids.count(pair1) == 1);
    REQUIRE (edge_ids.count(pair2) == 1);
    REQUIRE (edge_ids.count(pair3) == 1);
    REQUIRE (edge_ids.count(pair4) == 1);
    REQUIRE (edge_ids.count(pair5) == 1);
    REQUIRE (edge_ids.count(pair6) == 1);

    REQUIRE (edge_ids[pair1].to_ulong() == 3);
    REQUIRE (edge_ids[pair2].to_ulong() == 3);
    REQUIRE (edge_ids[pair3].to_ulong() == 1);
    REQUIRE (edge_ids[pair4].to_ulong() == 3);
    REQUIRE (edge_ids[pair5].to_ulong() == 2);
    REQUIRE (edge_ids[pair6].to_ulong() == 2);
  }
}


TEST_CASE ("extendGraph should extend the graph further")
{
  const char* alignment_file = (char*)"alignments/test_case5.fa";
  std::vector<std::string> ids;
  std::vector<VertexLabels> vertex_vector;
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;
  boost::unordered_set<TVertexDescriptor> free_nodes;
  free_nodes.insert(1);
  TVertexDescriptor new_begin_vertex;
  
  TVertexDescriptor begin_vertex;
  TGraph g = createGraph(alignment_file, vertex_vector, edge_ids, ids, begin_vertex);
  free_nodes.insert(begin_vertex);
  alignment_file = (char*)"../data/Alignments/test_case4_i2.fa";
  extendGraph(g, alignment_file, vertex_vector, edge_ids, new_begin_vertex, begin_vertex);

  SECTION ("graph should now have 2 more edges and 1 more vertix")
  {
    REQUIRE (numVertices(g) == 11);
    REQUIRE (numEdges(g) == 12);
  }

  SECTION ("there should be an edge between the labels")
  {
    String<bool> matrix;
    getAdjacencyMatrix(g, matrix);

    REQUIRE (length(matrix) == 11*11);
  }

  SECTION ("there shouldn't be any more edges in that matrix")
  {
    String<bool> matrix;
    getAdjacencyMatrix(g, matrix);

    unsigned sum = 0;
    for (unsigned i = 0 ; i < length(matrix) ; ++i)
    {
      if ( matrix[i] != 0 )
        ++sum;
    }
    REQUIRE (sum == 12);
  }
}
*/

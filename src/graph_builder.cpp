#include "graph_builder.hpp"

using namespace seqan;


inline void
addBitsetEdge (boost::dynamic_bitset<> &ref_bitset,
               const unsigned short &bit_id,
               const unsigned short &bit_n,
               const bool &non_empty_ref
              )
{
  if (non_empty_ref)
  {
    boost::dynamic_bitset<> new_bitset_edge(bit_n);
    new_bitset_edge[bit_id] = 1;
    ref_bitset = ref_bitset | new_bitset_edge;
  }
  else
  {
    boost::dynamic_bitset<> x(bit_n);
    x[bit_id] = 1;
    ref_bitset = x;
  }
}


void
addSequenceToGraph (TGraph & g,
                    unsigned short const & bit_id,
                    unsigned short const & bit_n,
                    CharString const & sequence,
                    std::map<VertexLabels, TVertexDescriptor> & vertex_label_map,
                    std::vector<VertexLabels> & vertex_vector,
                    boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids
                   )
{
  VertexLabels prev_vertex_label = {0, Dna('A')};
  TVertexDescriptor prev = vertex_label_map[prev_vertex_label];
  
  for (int pos = 0 ; ; ++pos)
  {
    if (pos >= (int) length(sequence))
    {
      VertexLabels final_vertex_label = {
        -1, // level
        Dna('A'), // Dna, initialized here but we should never use it!
      };

      TVertexDescriptor target_vertex = vertex_label_map[final_vertex_label];
      std::pair<TVertexDescriptor, TVertexDescriptor> my_pair (prev, target_vertex);
      addBitsetEdge(edge_ids[my_pair], bit_id, bit_n, edge_ids.count(my_pair));

      if (!findEdge(g, prev, target_vertex))
        addEdge(g, prev, target_vertex);

      break;
    }

    if (sequence[pos] == '-')
      continue;

    VertexLabels new_vertex_label =
    {
      pos+1,               // level
      Dna(sequence[pos]),  // dna
    };

    TVertexDescriptor target_vertex;

    if (vertex_label_map.count(new_vertex_label) == 1)
    {
      target_vertex = vertex_label_map[new_vertex_label];
      std::pair<TVertexDescriptor, TVertexDescriptor> my_pair (prev, target_vertex);
      addBitsetEdge(edge_ids[my_pair], bit_id, bit_n, edge_ids.count(my_pair));

      if (!findEdge(g, prev, target_vertex))
        addEdge(g, prev, target_vertex);
    }
    else
    {
      target_vertex = addVertex(g);
      std::pair<TVertexDescriptor, TVertexDescriptor> my_pair (prev, target_vertex);
      addBitsetEdge(edge_ids[my_pair], bit_id, bit_n, edge_ids.count(my_pair));
      addEdge(g, prev, target_vertex);
      vertex_label_map[new_vertex_label] = target_vertex;
      vertex_vector.push_back(new_vertex_label);
    }
    prev = target_vertex;
  }
}


void
addSequenceToGraph (TGraph & g,
                    CharString const & sequence,
                    std::map<VertexLabels, TVertexDescriptor> & vertex_label_map,
                    std::vector<VertexLabels> & vertex_vector
                   )
{
  // Used when we don't want to add to the edge_ids map.
  VertexLabels initial_vertex = { 0, Dna('A') };
  TVertexDescriptor prev = vertex_label_map[initial_vertex];
  
  for (int pos = 0 ; ; ++pos)
  {
    if (pos >= (int) length(sequence))
    {
      VertexLabels final_vertex_label = {
        -1, // level
        Dna('A'),
      };

      TVertexDescriptor target_vertex = vertex_label_map[final_vertex_label];

      if (!findEdge(g, prev, target_vertex))
      {
        addEdge(g, prev, target_vertex);
      }

      break;
    }

    if (sequence[pos] == '-')
      continue;

    VertexLabels new_vertex_label =
    {
      pos+1,               // level
      Dna(sequence[pos]),  // dna
    };

    if (vertex_label_map.count(new_vertex_label) == 1)
    {
      TVertexDescriptor target_vertex = vertex_label_map[new_vertex_label];

      if (!findEdge(g, prev, target_vertex))
      {
        addEdge(g, prev, target_vertex);
      }

      prev = target_vertex;
    }
    else
    {
      TVertexDescriptor target_vertex = addVertex(g);
      addEdge(g, prev, target_vertex);
      vertex_label_map[new_vertex_label] = target_vertex;
      prev = target_vertex;
      vertex_vector.push_back(new_vertex_label);
    }
  }
}


void
addInitialAndEndVertex (TGraph & g,
                        std::map<VertexLabels, TVertexDescriptor> & vertex_label_map,
                        std::vector<VertexLabels> & vertex_vector,
                        TVertexDescriptor & begin_vertex
                       )
{
  begin_vertex = addVertex(g);
  TVertexDescriptor end_vertex = addVertex(g);

  VertexLabels initial_vertex = 
  {
    0, // level
    Dna('A'), // Dna
  };

  VertexLabels final_vertex_label = {
    -1, // level
    Dna('A'), // Dna
  };

  vertex_label_map[initial_vertex] = begin_vertex;
  vertex_vector.push_back(initial_vertex);
  vertex_label_map[final_vertex_label] = end_vertex;
  vertex_vector.push_back(final_vertex_label);
}


void
addInitialVertex (TGraph &g,
                  std::map<VertexLabels, TVertexDescriptor> &vertex_label_map,
                  std::vector<VertexLabels> &vertex_vector,
                  TVertexDescriptor &begin_vertex,
                  TVertexDescriptor &end_vertex
                 )
{
  begin_vertex = addVertex(g);

  VertexLabels initial_vertex = 
  {
    0, // level
    Dna('A'), // Dna
  };

  VertexLabels final_vertex_label = {
    -1, // level
    Dna('A'), // Dna
  };

  vertex_label_map[initial_vertex] = begin_vertex;
  vertex_vector.push_back(initial_vertex);
  vertex_label_map[final_vertex_label] = end_vertex;
  // vertex_vector.push_back(final_vertex_label);

  end_vertex = begin_vertex;
}


TGraph
createGraph (const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > &edge_ids,
             std::vector<std::string> &ids,
             TVertexDescriptor &begin_vertex
            )
{
  TGraph graph;
  CharString id, sequence;
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  SeqFileIn seqFileIn(fastaFile);
  SeqFileIn countRecords(fastaFile);
  String<CharString> count_ids, count_seqs;
  readRecords(count_ids, count_seqs, countRecords);
  addInitialAndEndVertex(graph, vertex_label_map, vertex_vector, begin_vertex);

  unsigned id_read = 0;

  while (!atEnd(seqFileIn))
  {
    readRecord(id, sequence, seqFileIn);
    addSequenceToGraph(graph, id_read, length(count_ids), sequence, vertex_label_map, vertex_vector, edge_ids);
    std::string meta = toCString(id);
    ids.push_back(meta.substr(3));
    ++id_read;
  }

  return graph;
}


TGraph
createGraph (const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             std::vector<std::string> &ids,
             TVertexDescriptor &begin_vertex
            )
{
  // Used when we want to create a graph but no add anything to the edge_ids map.
  TGraph graph;
  CharString id, sequence;
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  SeqFileIn seqFileIn(fastaFile);

  addInitialAndEndVertex(graph, vertex_label_map, vertex_vector, begin_vertex);

  while (!atEnd(seqFileIn))
  {
    readRecord(id, sequence, seqFileIn);
    addSequenceToGraph(graph, sequence, vertex_label_map, vertex_vector);
    std::string meta = toCString(id);
    ids.push_back(meta.substr(3));
  }

  return graph;
}


void
extendGraph (TGraph &graph,
             const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > &edge_ids,
             TVertexDescriptor &begin_vertex,
             TVertexDescriptor &end_vertex
            )
{
  CharString id, sequence;
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  SeqFileIn seqFileIn(fastaFile);
  SeqFileIn countRecords(fastaFile);
  String<CharString> count_ids, count_seqs;
  readRecords(count_ids, count_seqs, countRecords);

  readRecord(id, sequence, seqFileIn);
  addLevelsToVertexVector(vertex_vector, length(sequence));
  addInitialVertex(graph, vertex_label_map, vertex_vector, begin_vertex, end_vertex);
  unsigned id_read = 0;

  while (!atEnd(seqFileIn))
  {
    addSequenceToGraph(graph, id_read, length(count_ids), sequence, vertex_label_map, vertex_vector, edge_ids);
    readRecord(id, sequence, seqFileIn);
    ++id_read;
  }
  
  addSequenceToGraph(graph, id_read, length(count_ids), sequence, vertex_label_map, vertex_vector, edge_ids);
}


void
extendGraph (TGraph &graph,
             const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             TVertexDescriptor &begin_vertex,
             TVertexDescriptor &end_vertex
            )
{
  // Used when graph is extended without added to the edge_ids map.
  CharString id, sequence;
  std::map<VertexLabels, TVertexDescriptor> vertex_label_map;
  SeqFileIn seqFileIn(fastaFile);

  readRecord(id, sequence, seqFileIn);
  addLevelsToVertexVector(vertex_vector, length(sequence));
  addInitialVertex(graph, vertex_label_map, vertex_vector, begin_vertex, end_vertex);

  while (!atEnd(seqFileIn))
  {
    addSequenceToGraph(graph, sequence, vertex_label_map, vertex_vector);
    readRecord(id, sequence, seqFileIn);
  }
  addSequenceToGraph(graph, sequence, vertex_label_map, vertex_vector);
}

#include <stdio.h>

#include "graph_align.hpp"
#include <cmath>
#include <set>
#include <climits>

using namespace seqan;


void
initializeExactScoreMatrixAndBacktracker(std::size_t const & seq_length,
                                         std::size_t const & number_of_nodes,
                                         std::vector<ExactBacktracker> & backtracker
                                        )
{
  backtracker.clear();
  backtracker.resize(number_of_nodes);

  for (TVertexDescriptor node = 0 ; node < number_of_nodes ; ++node)
  {
    std::vector< TVertexDescriptor > nodes;
    nodes.resize(seq_length, node);
    std::vector< bool > match;
    match.resize(seq_length, false);

    backtracker[node].nodes = nodes;
    backtracker[node].match = match;
  }
}


/*
 * getScoreVector 
 **/
  // With quality
bool
getScoreVector(ExactBacktracker const & previous_backtracker,
               ExactBacktracker & current_backtracker,
               DnaString const & sequence,
               boost::dynamic_bitset<> const & qual,
               Dna const & reference,
               TVertexDescriptor const & source_vertex
              )
{
  // I only care if there is a diagonal match when checking for exact matches only
  if (sequence[0] == reference || qual.test(0))
  {
    current_backtracker.match[0] = true;
    current_backtracker.nodes[0] = source_vertex;
  }
  
  std::size_t seq_length = length(sequence);

  for (unsigned pos = 1; pos < seq_length ; ++pos)
  {
    if (previous_backtracker.match[pos-1] && (sequence[pos] == reference || qual.test(pos)))
    {
      current_backtracker.match[pos] = true;
      current_backtracker.nodes[pos] = source_vertex;
    }
  }

  return current_backtracker.match.back();
}

// Without quality
bool
getScoreVector(ExactBacktracker const & previous_backtracker,
               ExactBacktracker & current_backtracker,
               DnaString const & sequence,
               Dna const & reference,
               TVertexDescriptor const & source_vertex
              )
{
  // I only care if there is a diagonal match when checking for exact matches only
  if (sequence[0] == reference)
  {
    current_backtracker.match[0] = true;
    current_backtracker.nodes[0] = source_vertex;
  }
  
  std::size_t seq_length = length(sequence);
  for (unsigned pos = 1; pos < seq_length ; ++pos)
  {
    if (previous_backtracker.match[pos-1] && sequence[pos] == reference)
    {
      current_backtracker.match[pos] = true;
      current_backtracker.nodes[pos] = source_vertex;
    }
  }

  return current_backtracker.match.back();
}

/*
 * alignToGraphExact aligns a sequence to the graph
 **/
// With quality
void
alignToGraphExact (DnaString const & sequence,
                   String<TVertexDescriptor const> const & order,
                   TGraph const & graph,
                   std::vector<TVertexDescriptor> & matching_vertices,
                   std::vector<VertexLabels> & vertex_vector,
                   std::vector<ExactBacktracker> & backtracker,
                   boost::unordered_set<TVertexDescriptor> const & free_nodes,
                   boost::dynamic_bitset<> const & qual
                  )
{
  for (Iterator<String<TVertexDescriptor const> const>::Type it = begin(order) ; it != end(order) ; ++it)
  {
    TVertexDescriptor const & source_vertex = *it;
    for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
    {
      TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);

      if (free_nodes.find(target_vertex) != free_nodes.end())
      {
        for (unsigned i = 0 ; i < length(sequence) ; ++i)
        {
          if (backtracker[source_vertex].match[i])
          {
            backtracker[target_vertex].match[i] = true;
            backtracker[target_vertex].nodes[i] = source_vertex;
          }
        }
      }
      else if (outDegree(graph, source_vertex) == 1 &&
               inDegree(graph, target_vertex) == 1
              )
      {
        Dna reference = vertex_vector[target_vertex].dna;
        if (getScoreVector(backtracker[source_vertex],
                           backtracker[target_vertex],
                           sequence,
                           qual,
                           reference,
                           source_vertex
                          )
           )
        {
          matching_vertices.push_back(target_vertex);
        }
      }
      else
      {
        Dna reference = vertex_vector[target_vertex].dna;
        if (getScoreVector(backtracker[source_vertex],
                           backtracker[target_vertex],
                           sequence,
                           reference,
                           source_vertex
                          )
           )
        {
          matching_vertices.push_back(target_vertex);
        }
      }
    }
  }
}


// Without quality
void
alignToGraphExact (DnaString const & sequence,
                   String<TVertexDescriptor const> const & order,
                   TGraph const & graph,
                   std::vector<TVertexDescriptor> & matching_vertices,
                   std::vector<VertexLabels> & vertex_vector,
                   std::vector<ExactBacktracker> & backtracker,
                   boost::unordered_set<TVertexDescriptor> const & free_nodes
                  )
{
  for (Iterator<String<TVertexDescriptor const> const>::Type it = begin(order) ; it != end(order) ; ++it)
  {
    TVertexDescriptor const & source_vertex = *it;
    for (Iterator<TGraph, OutEdgeIterator>::Type out_edge_iterator (graph, source_vertex) ; !atEnd(out_edge_iterator) ; ++out_edge_iterator)
    {
      TVertexDescriptor const & target_vertex = targetVertex(out_edge_iterator);

      if (free_nodes.find(target_vertex) != free_nodes.end())
      {
        for (unsigned i = 0 ; i < length(sequence) ; ++i)
        {
          if (backtracker[source_vertex].match[i])
          {
            backtracker[target_vertex].match[i] = true;
            backtracker[target_vertex].nodes[i] = source_vertex;
          }
        }
      }
      else
      {
        Dna reference = vertex_vector[target_vertex].dna;
        if (getScoreVector(backtracker[source_vertex],
                           backtracker[target_vertex],
                           sequence,
                           reference,
                           source_vertex
                          )
           )
        {
          matching_vertices.push_back(target_vertex);
        }
      }
    }
  }
}


boost::dynamic_bitset<>
backTrackAndCount (std::vector<std::string> const & ids,
                   std::vector<ExactBacktracker> const & backtracker,
                   TVertexDescriptor old_node_id,
                   boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids
                  )
{
  // #define __backTrackAndCount_VERBOSE__
  boost::unordered_map<std::string, long> current_map;
  boost::dynamic_bitset<> perfect_match(length(ids));
  perfect_match = ~perfect_match;

  for (unsigned i = 0 ; i < ids.size() ; ++i)
  {
    current_map[ids[i]] = 0;
  }

  TVertexDescriptor new_node_id;
  int row = backtracker[old_node_id].nodes.size()-1;

  #ifdef __backTrackAndCount_VERBOSE__
  std::cout << "perfect_match = " << perfect_match << std::endl;
  #endif

  while (row >= 0)
  {
    new_node_id = backtracker[old_node_id].nodes[row];
    #ifdef __backTrackAndCount_VERBOSE__
    std::cout << "old_node_id = " << old_node_id << " -> " << "new_node_id = " << new_node_id << " with row = " << row << std::endl;
    std::cout << "backtracker[" << old_node_id << "].match[" << row << "] = " << backtracker[old_node_id].match[row] << std::endl;
    #endif

    if (backtracker[old_node_id].match[row])
    {
      // std::string key_string = boost::lexical_cast<std::string>(new_node_id).append(">").append(boost::lexical_cast<std::string>(old_node_id));

      std::pair<TVertexDescriptor, TVertexDescriptor> key_pair(new_node_id, old_node_id);

      #ifdef __backTrackAndCount_VERBOSE__
      std::cout << "We have a match!" << std::endl;
      #endif

      old_node_id = new_node_id;

      --row;
      if (edge_ids.count(key_pair) == 0 || row == -1)
      {
        // We don't want to count edges that go to row -1 cause we could've picked any edge for that! This is a neglatable bias when graph is big.
        #ifdef __backTrackAndCount_VERBOSE__
        std::cout << "Edge is on a intron." << std::endl;
        #endif

        continue;
      }

      perfect_match = perfect_match & edge_ids[key_pair];
    }
    else
    {
      ++row;
      #ifdef __backTrackAndCount_VERBOSE__
      std::cout << "Walked over a free node." << std::endl;
      #endif
    }
  }

  #ifdef __backTrackAndCount_VERBOSE__
  std::cout << "perfect_match = " << perfect_match << std::endl;
  #endif

  return perfect_match;
}

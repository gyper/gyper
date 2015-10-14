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


boost::dynamic_bitset<>
align_kmer_to_graph (String<Dna> const & sequence,
                     unsigned const & id_numbers,
                     TKmerMap & kmer_map,
                     unsigned const & best_kmer_index,
                     int const & kmer_size
                    )
{
  String<Dna> seq_best_kmer(sequence);
  erase(seq_best_kmer, 0, best_kmer_index);
  resize(seq_best_kmer, kmer_size);
  std::vector< KmerLabels > matches;
  boost::dynamic_bitset<> id_bits(id_numbers);

  if (kmer_map.count(seq_best_kmer) == 0)
  {
    return id_bits;
  }
  else
  {
    matches = kmer_map[seq_best_kmer];
  }

  int const & increment_size = kmer_size - 1;

  // TODO: Make this loop simpler
  for (unsigned k = best_kmer_index + increment_size ; k < length(sequence) - increment_size ; k += increment_size)
  {
    // std::cout << "Sequence: " << sequence << std::endl;
    // std::cout << "matches.size() = " << matches.size() << std::endl;
    String<Dna> seq_center_kmer(sequence);
    erase(seq_center_kmer, 0, k);
    resize(seq_center_kmer, kmer_size);

    if (kmer_map.count(seq_center_kmer) == 1)
    {
      // std::cout << "Match for kmer " << seq_center_kmer << std::endl;
      unsigned new_matches = 0;

      for (unsigned pos = 0; pos < matches.size() - new_matches ; ++pos)
      {
        KmerLabels original_matches = matches[pos];
        // std::cout << "pos = " << pos << " is searching for " << matches[pos].end_vertex << std::endl;
        bool matched = false;

        for (auto current_matches_it = kmer_map[seq_center_kmer].begin() ; current_matches_it != kmer_map[seq_center_kmer].end() ; ++current_matches_it)
        {
          // std::cout << "Considering kmer " << seq_center_kmer << ": " << current_matches_it->start_vertex << " " << current_matches_it->end_vertex << " " << current_matches_it->id_bits << std::endl;

          if (original_matches.end_vertex == current_matches_it->start_vertex)
          {
            if (matched)
            {
              // std::cout << "Another match found! " << current_matches_it->start_vertex << " " << current_matches_it->end_vertex << " " << current_matches_it->id_bits << std::endl;
              KmerLabels new_label = {
                original_matches.start_vertex,
                current_matches_it->end_vertex,
                original_matches.id_bits & current_matches_it->id_bits
              };

              matches.push_back(new_label);
              ++new_matches;
            }
            else
            {
              // std::cout << "Match found! " << current_matches_it->start_vertex << " " << current_matches_it->end_vertex << " " << current_matches_it->id_bits << std::endl;
              matched = true;
              matches[pos].end_vertex = current_matches_it->end_vertex;
              matches[pos].id_bits &= current_matches_it->id_bits;
            }
          }
        }

        if (!matched)
        {
          // std::cout << "Erased. " << matches[pos].start_vertex << " " << matches[pos].end_vertex << " " << matches[pos].id_bits << std::endl;
          matches.erase(matches.begin() + pos);
          --pos;
        }
      }
 
      // matches = new_matches;
      // std::cout << "Kmer: " << seq_center_kmer << " (" << matches.size() << ") " << std::endl;
    }
    else
    {
      return id_bits;
    }

    if (matches.size() == 0)
    {
      return id_bits;
    }
  }

  // std::cout << "Sequence I have is " << sequence << " with " << matches.size() << " kmers." << std::endl;

  for (auto matches_it = matches.begin() ; matches_it != matches.end() ; ++matches_it)
  {
    id_bits |= matches_it->id_bits;
  }

  return id_bits;
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

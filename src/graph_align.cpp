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


bool
find_forward_matches(std::vector< KmerLabels > & matches,
                     String<Dna> const & kmer,
                     TKmerMap & kmer_map
                    )
{
  std::vector< KmerLabels > original_matches(matches);
  unsigned new_matches = 0;

  for (unsigned pos = 0; pos < matches.size() - new_matches ; ++pos)
  {
    KmerLabels original_matches_pos = matches[pos];
    bool matched = false;

    for (auto current_matches_it = kmer_map[kmer].begin() ; current_matches_it != kmer_map[kmer].end() ; ++current_matches_it)
    {
      if (original_matches_pos.end_vertex == current_matches_it->start_vertex)
      {
        if (matched)
        {
          // std::cout << "Another match found! " << current_matches_it->start_vertex << " " << current_matches_it->end_vertex << " " << current_matches_it->id_bits << std::endl;
          KmerLabels new_label = {
            original_matches_pos.start_vertex,
            current_matches_it->end_vertex,
            original_matches_pos.id_bits & current_matches_it->id_bits
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
  
  if (matches.size() == 0)
  {
    matches = original_matches;
    return false;
  }

  return true;
}


bool
find_backward_matches(std::vector< KmerLabels > & matches,
                      String<Dna> const & kmer,
                      TKmerMap & kmer_map
                     )
{
  std::vector< KmerLabels > original_matches(matches);
  unsigned new_matches = 0;

  for (unsigned pos = 0; pos < matches.size() - new_matches ; ++pos)
  {
    KmerLabels original_matches_pos = matches[pos];
    bool matched = false;

    for (auto current_matches_it = kmer_map[kmer].begin() ; current_matches_it != kmer_map[kmer].end() ; ++current_matches_it)
    {
      if (original_matches_pos.start_vertex == current_matches_it->end_vertex)
      {
        if (matched)
        {
          // std::cout << "Another match found! " << current_matches_it->start_vertex << " " << current_matches_it->end_vertex << " " << current_matches_it->id_bits << std::endl;
          KmerLabels new_label = {
            current_matches_it->start_vertex,
            original_matches_pos.end_vertex,
            original_matches_pos.id_bits & current_matches_it->id_bits
          };

          matches.push_back(new_label);
          ++new_matches;
        }
        else
        {
          // std::cout << "Match found! " << current_matches_it->start_vertex << " " << current_matches_it->end_vertex << " " << current_matches_it->id_bits << std::endl;
          matched = true;
          matches[pos].start_vertex = current_matches_it->start_vertex;
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
  
  if (matches.size() == 0)
  {
    matches = original_matches;
    return false;
  }

  return true;
}


boost::dynamic_bitset<>
align_kmer_to_graph (String<Dna> const & sequence,
                     unsigned const & id_numbers,
                     TKmerMap & kmer_map,
                     std::vector<VertexLabels> & vertex_vector,
                     unsigned const & best_kmer_index,
                     int const & kmer_size,
                     int min_kmers
                    )
{
  String<Dna> seq_best_kmer(sequence);
  erase(seq_best_kmer, 0, best_kmer_index);

  // Debugging
  // if (length(seq_best_kmer) < static_cast<std::size_t>(kmer_size))
  // {
  //   std::cerr << "Best kmer is shorter than k! " << seq_best_kmer << std::endl;
  // }

  resize(seq_best_kmer, kmer_size);
  std::vector< KmerLabels > matches;
  boost::dynamic_bitset<> id_bits(id_numbers);

  if (kmer_map.count(seq_best_kmer) == 0)
  {
    // If we cant find the best_kmer, just give up
    // std::cout << "Didn't find " << seq_best_kmer << std::endl;
    return id_bits;
  }
  else
  {
    matches = kmer_map[seq_best_kmer];
  }

  --min_kmers;

  // Debugging
  // boost::dynamic_bitset<> id_bits_debug(id_numbers);
  // for (auto matches_it = matches.begin() ; matches_it != matches.end() ; ++matches_it)
  // {
  //   id_bits_debug |= matches_it->id_bits;
  // }
  // std::cout << "id_bits, best kmer = " << id_bits_debug << " " << seq_best_kmer << std::endl;

  int const & increment_size = kmer_size - 1;

  bool check_forward_potential = true;

  // Forward loop
  for (int k = best_kmer_index + increment_size ; k < static_cast<int>(length(sequence)) - increment_size ; k += increment_size)
  {
    String<Dna> kmer(sequence);
    erase(kmer, 0, k);
    resize(kmer, kmer_size);

    // std::cout << kmer << std::endl;

    if (kmer_map.count(kmer) == 0 || !find_forward_matches(matches, kmer, kmer_map))
    {
      // std::cout << "Mismatch forward." << std::endl;
      check_forward_potential = false;
      break;
    }

    --min_kmers;
  }

  boost::dynamic_bitset<> potential_forward(id_numbers);

  if (check_forward_potential)
  {
    String<Dna> seq_last_kmer(sequence);
    erase(seq_last_kmer, 0, length(sequence) - kmer_size);
    
    if (kmer_map.count(seq_last_kmer) == 1)
    {
      for (auto current_matches_it = kmer_map[seq_last_kmer].begin() ; current_matches_it != kmer_map[seq_last_kmer].end() ; ++current_matches_it)
      {
        for (auto original_matches_it = matches.begin() ; original_matches_it != matches.end() ; ++original_matches_it)
        {
          if (abs(vertex_vector[current_matches_it->start_vertex].level - vertex_vector[original_matches_it->start_vertex].level) < kmer_size * 2)
          {
            potential_forward |= current_matches_it->id_bits;
            continue;
          }
        }
      }
    }
  }


  bool check_backward_potential = true;

  // Backward loop
  for (int k = best_kmer_index - increment_size ; k >= 0 ; k -= increment_size)
  {
    String<Dna> kmer(sequence);
    erase(kmer, 0, k);
    resize(kmer, kmer_size);

    // std::cout << kmer << std::endl;

    if (kmer_map.count(kmer) == 0 || !find_backward_matches(matches, kmer, kmer_map))
    {
      // std::cout << "Mismatch backward." << std::endl;
      check_backward_potential = false;
      break;
    }

    --min_kmers;
  }

  boost::dynamic_bitset<> potential_backward(id_numbers);

  // Check if the last kmer could potentially give any information
  if (check_backward_potential)
  {
    // std::cout << "Checking backward potentially" << std::endl;
    String<Dna> seq_first_kmer(sequence);
    resize(seq_first_kmer, kmer_size);

    if (kmer_map.count(seq_first_kmer) == 1)
    {
      for (auto current_matches_it = kmer_map[seq_first_kmer].begin() ; current_matches_it != kmer_map[seq_first_kmer].end() ; ++current_matches_it)
      {
        for (auto original_matches_it = matches.begin() ; original_matches_it != matches.end() ; ++original_matches_it)
        {
          if (abs(vertex_vector[current_matches_it->start_vertex].level - vertex_vector[original_matches_it->start_vertex].level) < kmer_size * 2)
          {
            potential_backward |= current_matches_it->id_bits;
            continue;
          }
        }
      }
    }
  }

  if (min_kmers > 0)
  {
    return id_bits;
  }
  
  for (auto matches_it = matches.begin() ; matches_it != matches.end() ; ++matches_it)
  {
    id_bits |= matches_it->id_bits;
  }
  

  if (potential_forward.none())
  {
    potential_forward.flip();
  }
  else
  {
    --min_kmers;
  }

  if (potential_backward.none())
  {
    potential_backward.flip();
  }
  else
  {
    --min_kmers;
  }

  return id_bits & potential_forward & potential_backward;
  // return id_bits;
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

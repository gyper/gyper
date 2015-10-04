#ifndef __GRAPH_ALIGN_HPP_INCLUDED__
#define __GRAPH_ALIGN_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include "graph_builder.hpp"

#include <boost/unordered_set.hpp>

struct ExactBacktracker {
  std::vector<TVertexDescriptor> nodes;
  std::vector<bool> match; // True if diagonal match, false otherwise
};

#include "graph_io.hpp"
// Needs to be here cause it's dependent on Backtracker and ExactBacktracker

void
initializeExactScoreMatrixAndBacktracker(std::size_t const & seq_length,
                                         std::size_t const & number_of_nodes,
                                         std::vector<ExactBacktracker> & backtracker
                                        );

bool
getScoreVector(ExactBacktracker const & previous_backtracker,
               ExactBacktracker & current_backtracker,
               DnaString const & sequence,
               boost::dynamic_bitset<> const & qual,
               Dna const & reference,
               TVertexDescriptor const & source_vertex
              );

bool
getScoreVector(ExactBacktracker const & previous_backtracker,
               ExactBacktracker & current_backtracker,
               DnaString const & sequence,
               Dna const & reference,
               TVertexDescriptor const & source_vertex
              );

void
alignToGraphExact (DnaString const & sequence,
                   String<TVertexDescriptor const> const & order,
                   TGraph const & graph,
                   std::vector<TVertexDescriptor> & matching_vertices,
                   std::vector<VertexLabels> & vertex_vector,
                   std::vector<ExactBacktracker> & backtracker,
                   boost::unordered_set<TVertexDescriptor> const & free_nodes,
                   boost::dynamic_bitset<> const & qual
                  );

void
alignToGraphExact (DnaString const & sequence,
                   String<TVertexDescriptor const> const & order,
                   TGraph const & graph,
                   std::vector<TVertexDescriptor> & matching_vertices,
                   std::vector<VertexLabels> & vertex_vector,
                   std::vector<ExactBacktracker> & backtracker,
                   boost::unordered_set<TVertexDescriptor> const & free_nodes
                  );

boost::dynamic_bitset<>
backTrackAndCount (std::vector<std::string> const & ids,
                   std::vector<ExactBacktracker> const & backtracker,
                   TVertexDescriptor old_node_id,
                   boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids
                  );

#endif

#ifndef __GRAPH_HPP_INCLUDED__
#define __GRAPH_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include <stdio.h>
#include <ctime>

#include "graph_align.hpp"
#include "graph_builder.hpp"
#include "graph_io.hpp"

#include <boost/unordered/unordered_set.hpp>

#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include <seqan/seq_io.h>


using namespace seqan;


void
align_sequence (DnaString & my_sequence,
                boost::dynamic_bitset<> & qual,
                TGraph const & graph,
                std::vector<VertexLabels> & vertex_vector,
                String<TVertexDescriptor> & order,
                std::vector<ExactBacktracker> & backtracker,
                std::vector<ExactBacktracker> & reverse_backtracker,
                boost::unordered_set<TVertexDescriptor> const & free_nodes,
                std::vector<TVertexDescriptor> & matching_vertices,
                std::vector<TVertexDescriptor> & reverse_matching_vertices
               );

CharString myExtractTagValue(String<char> &tags);

void
createDqa1Graph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createDqb1Graph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createDrb1Graph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createHlaaGraph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createHlabGraph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createHlacGraph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

#endif

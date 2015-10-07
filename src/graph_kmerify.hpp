#ifndef __GRAPH_KMERIFY_HPP_INCLUDED__
#define __GRAPH_KMERIFY_HPP_INCLUDED__
#include <string>
#include <sstream>
#include <iostream>

#include "graph_builder.hpp"
#include "graph_align.hpp"

#include <boost/unordered_set.hpp>

#define K_SIZE 15

boost::unordered_map< std::string, std::vector<TVertexDescriptor> >
kmerifyGraph(String<TVertexDescriptor const> const & order,
             TGraph const & graph,
             std::vector<VertexLabels> & vertex_vector,
             boost::unordered_set<TVertexDescriptor> const & free_nodes
            );

#endif

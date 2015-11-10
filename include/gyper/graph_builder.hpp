#ifndef __GRAPH_BUILDER_HPP_INCLUDED__
#define __GRAPH_BUILDER_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include <gyper/graph.hpp>


inline bool
hasEdge (seqan::String<bool> adjacency_matrix, unsigned id_from, unsigned id_to, unsigned n)
{
  return adjacency_matrix[id_to + n * id_from];
}


inline bool
hasEdge (seqan::String<bool> adjacency_matrix, unsigned id_from, unsigned id_to)
{
  return hasEdge(adjacency_matrix, id_from, id_to, std::sqrt(length(adjacency_matrix)) );
}


inline void
addLevelsToVertexVector (std::vector<VertexLabels> & vertex_vector,
                         std::size_t const & level
                        )
{
  for (std::vector<VertexLabels>::iterator it = vertex_vector.begin() ; it != vertex_vector.end() ; ++it)
  {
    it->level += level;
  }
}


void
addExonToGraph (TGraph & g,
                unsigned short const & bit_id,
                unsigned short const & bit_n,
                seqan::CharString const & sequence,
                std::map<VertexLabels, TVertexDescriptor> & vertex_label_map,
                std::vector<VertexLabels> & vertex_vector,
                boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids
               );

void
addSequenceToGraph (TGraph & g,
                    seqan::CharString const & sequence,
                    std::map<VertexLabels, TVertexDescriptor> & vertex_label_map,
                    std::vector<VertexLabels> & vertex_vector
                   );

void
addInitialAndEndVertex (TGraph & g,
                        std::map<VertexLabels, TVertexDescriptor> & vertex_label_map,
                        std::vector<VertexLabels> & vertex_vector,
                        TVertexDescriptor & begin_vertex
                       );

void
addInitialVertex (TGraph &g,
                  std::map<VertexLabels, TVertexDescriptor> &vertex_label_map,
                  std::vector<VertexLabels> &vertex_vector,
                  TVertexDescriptor &begin_vertex,
                  TVertexDescriptor &end_vertex
                 );

TGraph
createGraph (const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > &edge_ids,
             std::vector<std::string> &ids,
             TVertexDescriptor &begin_vertex
            );

TGraph
createGraph (const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             std::vector<std::string> &ids,
             TVertexDescriptor &begin_vertex
            );

void
extendGraph (TGraph &graph,
             const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > &edge_ids,
             TVertexDescriptor &begin_vertex,
             TVertexDescriptor &end_vertex
            );

void
extendGraph (TGraph &graph,
             const char* fastaFile,
             std::vector<VertexLabels> &vertex_vector,
             TVertexDescriptor &begin_vertex,
             TVertexDescriptor &end_vertex
            );

#endif

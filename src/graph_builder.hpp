#ifndef __GRAPH_BUILDER_HPP_INCLUDED__
#define __GRAPH_BUILDER_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include <cstddef>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/seq_io.h>

#include <boost/unordered/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/dynamic_bitset.hpp>


using namespace seqan;


struct VertexLabels {
  int level;
  Dna dna;
};

// For graph
typedef Graph<Directed<void, WithSourceId> > TGraph;
typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
typedef Size<TGraph>::Type TSize;

// For property maps
typedef String<unsigned> TLevel;
typedef String<TLevel> TProperties;

typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

struct GraphWithLabels {
  TGraph graph;
  std::vector<VertexLabels> vertex_labels;
};


inline bool
operator==(const VertexLabels &lhs, const VertexLabels &rhs)
{
  return lhs.level == rhs.level && lhs.dna == rhs.dna;
}


inline bool
operator<(const VertexLabels &lhs, const VertexLabels &rhs)
{
  return lhs.level < rhs.level || (lhs.level == rhs.level && lhs.dna < rhs.dna);
}


inline std::size_t
hash_value(const VertexLabels &a)
{
  std::size_t seed = 0;
  boost::hash_combine(seed, a.level);
  boost::hash_combine(seed, ordValue(a.dna));
  return seed;
}


inline bool
hasEdge (String<bool> adjacency_matrix, unsigned id_from, unsigned id_to, unsigned n)
{
  return adjacency_matrix[id_to + n * id_from];
}


inline bool
hasEdge (String<bool> adjacency_matrix, unsigned id_from, unsigned id_to)
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
addSequenceToGraph (TGraph & g,
                    unsigned short const & bit_id,
                    unsigned short const & bit_n,
                    CharString const & sequence,
                    std::map<VertexLabels, TVertexDescriptor> & vertex_label_map,
                    std::vector<VertexLabels> & vertex_vector,
                    boost::unordered_map<std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids
                   );

void
addSequenceToGraph (TGraph & g,
                    CharString const & sequence,
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

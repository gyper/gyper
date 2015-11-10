#ifndef __GRAPH_CREATION_HPP_INCLUDED__
#define __GRAPH_CREATION_HPP_INCLUDED__

#include <seqan/basic.h>
#include <seqan/sequence.h>

class VertexLabel
{
 public:
  unsigned order;
  seqan::String<seqan::Dna> dna;

  VertexLabel();
  VertexLabel(unsigned order, seqan::String<seqan::Dna> dna);
};

inline bool
operator==(const VertexLabel &lhs, const VertexLabel &rhs)
{
  return lhs.order == rhs.order && lhs.dna == rhs.dna;
}


inline bool
operator<(const VertexLabel &lhs, const VertexLabel &rhs)
{
  return lhs.order < rhs.order || (lhs.order == rhs.order && lhs.dna < rhs.dna);
}

#endif

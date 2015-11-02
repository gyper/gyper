#ifndef __GRAPH_CREATION_HPP_INCLUDED__
#define __GRAPH_CREATION_HPP_INCLUDED__

#include <seqan/basic.h>
#include <seqan/sequence.h>

class VertexLabel
{
 public:
  unsigned level;
  seqan::String<seqan::Dna> dna;

  VertexLabel();
  VertexLabel(unsigned level, seqan::String<seqan::Dna> dna);
};

inline bool
operator==(const VertexLabel &lhs, const VertexLabel &rhs)
{
  return lhs.level == rhs.level && lhs.dna == rhs.dna;
}


inline bool
operator<(const VertexLabel &lhs, const VertexLabel &rhs)
{
  return lhs.level < rhs.level || (lhs.level == rhs.level && lhs.dna < rhs.dna);
}

#endif

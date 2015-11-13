#ifndef __GRAPH_CREATION_HPP_INCLUDED__
#define __GRAPH_CREATION_HPP_INCLUDED__

#include <seqan/basic.h>
#include <seqan/sequence.h>

namespace gyper
{

class VertexLabel
{
 public:
  uint32_t order;
  seqan::String<seqan::Dna> dna;
  unsigned variant_id;

  VertexLabel()
    : order(0), dna(""), variant_id(0) {}

  VertexLabel(uint32_t const & order, seqan::String<seqan::Dna> const & dna)
    : order(order), dna(dna), variant_id(0) {}

  VertexLabel(uint32_t const & order, seqan::String<seqan::Dna> const & dna, unsigned const & variant_id)
    : order(order), dna(dna), variant_id(variant_id) {}
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

} // namespace gyper

#endif

#include <gyper/graph_creation.hpp>

VertexLabel::VertexLabel()
{
  order = 0;
}

VertexLabel::VertexLabel(unsigned order, seqan::String<seqan::Dna> dna)
{
  VertexLabel::order = order;
  VertexLabel::dna = dna;
}

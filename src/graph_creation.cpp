#include "graph_creation.hpp"

VertexLabel::VertexLabel()
{
  level = 0;
}

VertexLabel::VertexLabel(unsigned level, seqan::String<seqan::Dna> dna)
{
  VertexLabel::level = level;
  VertexLabel::dna = dna;
}
#include "partial_order_graph.hpp"

void POGraph::set_values (int x, int y) {
  width = x;
  height = y;
}

int POGraph::area ()
{
  return POGraph::width*POGraph::height;
}


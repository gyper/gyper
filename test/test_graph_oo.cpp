#include <stdio.h>
#include <climits>

#include "../src/partial_order_graph.hpp"
#include "../src/graph.hpp"
#include "catch.hpp"
#include "constants.hpp"

TEST_CASE("POGraph object should have a constructor")
{
  callOptions CO = callOptions(); // Default call options


  POGraph test_pog = POGraph();
  REQUIRE(numEdges(test_pog.graph) > 0);
}
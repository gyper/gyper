#include <stdio.h>
#include <climits>

#include "../src/partial_order_graph.hpp"
#include "../src/graph.hpp"
#include "catch.hpp"
#include "constants.hpp"

TEST_CASE("Gyper object should have a working constructor")
{
  SECTION("Gyper object should have a default constructor ")
  {
    Gyper gyper = Gyper();
    REQUIRE(numEdges(gyper.graph) > 0);
    REQUIRE(gyper.CO.gene == "DQA1");
  }

  SECTION("Gyper object should have constructor which takes argument options as input")
  {
    callOptions CO = callOptions(); // Default call options
    CO.gene = "DQB1";
    Gyper gyper = Gyper(CO);
    REQUIRE(numEdges(gyper.graph) > 0);
    REQUIRE(gyper.CO.gene == "DQB1");
  }
}
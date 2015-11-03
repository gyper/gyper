#include <stdio.h>
#include <climits>

#include "../src/partial_order_graph.hpp"
// #include "../src/graph.hpp"
#include "catch.hpp"
#include "constants.hpp"


// TEST_CASE("Gyper object should have a working constructor")
// {
//   SECTION("Gyper object should have a default constructor ")
//   {
//     Gyper gyper = Gyper();
//     REQUIRE(numEdges(gyper.graph) == 0);
    
//     // Create a graph for DQB1
//     gyper.create_HLA_graph();
//     REQUIRE(numEdges(gyper.graph) > 0);
//     REQUIRE(gyper.CO.gene == "DQA1");
//   }

//   SECTION("Gyper object should have constructor which takes argument options as input")
//   {
//     Options CO = Options();
//     CO.gene = "DQB1";
//     Gyper gyper = Gyper(CO);
//     REQUIRE(numEdges(gyper.graph) == 0);

//     // Create a graph for DQB1
//     gyper.create_HLA_graph();
//     REQUIRE(numEdges(gyper.graph) > 0);
//     REQUIRE(gyper.CO.gene == "DQB1");
//   }
// }

TEST_CASE("Create a graph from a reference FASTA file")
{
  // Get the base path to test reference genomes
  std::stringstream base_path;
  base_path << gyper_SOURCE_DIRECTORY << "/test/reference/reference.fa";
  std::string base_path_str = base_path.str();

  SECTION("using all sequences in the FASTA")
  {
    Options CO = Options();

    seqan::String<char> reference_file_name_1 = base_path_str;
    FastaHackAPI reference_fasta = FastaHackAPI(reference_file_name_1);
    Gyper gyper = Gyper(CO, reference_fasta);
    std::string region = "chr1";
    REQUIRE(numEdges(gyper.graph) == 0);
    gyper.create_reference_graph(region);
    std::cout << "numEdges(gyper.graph) = " << numEdges(gyper.graph) << std::endl;
    REQUIRE(numEdges(gyper.graph) > 0);
    // std::cout << "Graph: " << gyper.graph << std::endl;
  }
}
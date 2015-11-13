#include <catch.hpp>

#include <stdio.h>
#include <climits>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>

#include <gyper/constructor.hpp>
#include <gyper/options.hpp>
#include <gyper/constants.hpp>
#include <gyper/indexer.hpp>


gyper::Indexer*
BEFORE(const char * region)
{
  std::stringstream vcf_path;
  vcf_path << gyper_SOURCE_DIRECTORY << "/test/reference/indexz.vcf.gz";
  std::string vcf_path_str = vcf_path.str();

  std::stringstream reference_path;
  reference_path << gyper_SOURCE_DIRECTORY << "/test/reference/index_test.fa";
  std::string reference_path_str = reference_path.str();

  gyper::Options* CO = new gyper::Options();
  gyper::Constructor* constructor = new gyper::Constructor(*CO);
  delete CO;

  constructor->construct_graph(reference_path_str.c_str(), vcf_path_str.c_str(), region);

  gyper::Indexer* indexer = new gyper::Indexer(constructor);
  // indexer->index_graph();
  // std::cout << "indexer->k = " << indexer->k << std::endl;
  // std::cout << "indexer->current_mer = " << indexer->current_mer << std::endl;
  return indexer;
}

// TEST_CASE("Constructor")
// {
//   gyper::Indexer* indexer = BEFORE();
//   REQUIRE(indexer->k == 32);
// }

TEST_CASE("FIRST STEPS")
{
  gyper::Indexer* indexer = BEFORE("chr1");
  // std::cout << "Graph:" << indexer->c->graph << std::endl;
  indexer->index_graph();
}

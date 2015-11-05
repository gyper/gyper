#include <stdio.h>
#include <climits>

#include "../src/partial_order_graph.hpp"
#include "../src/gyper_options.hpp"
// #include "../src/graph.hpp"
#include "catch.hpp"
#include "constants.hpp"

#include <seqan/stream.h>


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

TEST_CASE("FASTA I/O")
{
  // Get the base path to test reference genomes
  std::stringstream base_path;
  base_path << gyper_SOURCE_DIRECTORY << "/test/reference/reference.fa";
  std::string base_path_str = base_path.str();

  SECTION("Extract subsequences from a FASTA file.")
  {
    Options* CO = new Options();
    Gyper* gyper = new Gyper(*CO);
    delete CO;
    gyper->open_fasta(base_path_str.c_str());
    unsigned fasta_index_id = gyper->get_fasta_index_id("chr1");
    REQUIRE(fasta_index_id == 0);
    fasta_index_id = gyper->get_fasta_index_id("chr2");
    REQUIRE(fasta_index_id == 1);
    seqan::String<seqan::Dna5> sequence_infix;
    seqan::readRegion(sequence_infix, gyper->fasta_index, fasta_index_id, 0, 1);
    REQUIRE(sequence_infix == "N");
    seqan::readRegion(sequence_infix, gyper->fasta_index, fasta_index_id, 0, 100);
    REQUIRE(sequence_infix == "NNAGGTTTCCCCAGGTTTCCCC");
    seqan::readRegion(sequence_infix, gyper->fasta_index, fasta_index_id, 1, 0);
    REQUIRE(sequence_infix == "");

    delete gyper;
  }

  // SECTION("using all sequences in the FASTA.")
  // {
  //   Options CO = Options();

  //   seqan::String<char> reference_file_name_1 = base_path_str;
  //   FastaHackAPI reference_fasta = FastaHackAPI(reference_file_name_1);
  //   Gyper gyper = Gyper(CO, reference_fasta);
  //   std::string region = "chr1";
  //   REQUIRE(numEdges(gyper.graph) == 0);
  //   gyper.create_reference_graph(region);
  //   std::cout << "numEdges(gyper.graph) = " << numEdges(gyper.graph) << std::endl;
  //   REQUIRE(numEdges(gyper.graph) > 0);
  //   // std::cout << "Graph: " << gyper.graph << std::endl;
  // }
}

TEST_CASE("VCF I/O")
{
  // Get the base path to test reference genomes
  SECTION("Extract subsequences from a VCF file.")
  {
    std::stringstream base_path;
    base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example.vcf";
    std::string base_path_str = base_path.str();

    Options CO = Options();
    Gyper* gyper = new Gyper(CO);
    gyper->open_vcf(base_path_str.c_str());

    REQUIRE(gyper->read_vcf_record() == 0);
    REQUIRE(gyper->vcf_record.rID == 0);
    REQUIRE(gyper->vcf_record.beginPos == 1);
    REQUIRE(gyper->vcf_record.id == "rs6054257");
    REQUIRE(gyper->vcf_record.ref == "G");
    REQUIRE(gyper->vcf_record.alt == "A");
    REQUIRE(gyper->vcf_record.qual == 0.0);
    REQUIRE(gyper->vcf_record.filter == ".");
    REQUIRE(gyper->vcf_record.info == ".");
    REQUIRE(gyper->vcf_record.format == "");
    REQUIRE(length(gyper->vcf_record.genotypeInfos) == 0);
    seqan::String<seqan::Dna> ref_dna = "G";
    REQUIRE(static_cast<seqan::String<seqan::Dna>>(gyper->vcf_record.ref) == ref_dna);

    REQUIRE(gyper->read_vcf_record() == 0);
    REQUIRE(gyper->vcf_record.rID == 0);
    REQUIRE(gyper->vcf_record.id == "ins504320");
    seqan::String<seqan::Dna> alt_dna = "GT";
    REQUIRE(static_cast<seqan::String<seqan::Dna>>(gyper->vcf_record.alt) == alt_dna);

    REQUIRE(gyper->read_vcf_record() == 0);
    REQUIRE(gyper->vcf_record.rID == 0);
    REQUIRE(gyper->vcf_record.id == "del504320");

    REQUIRE(gyper->read_vcf_record() == 0);
    REQUIRE(gyper->vcf_record.rID == 0);
    REQUIRE(gyper->vcf_record.id == "msat14123");

    REQUIRE(gyper->read_vcf_record() == 0);
    REQUIRE(gyper->vcf_record.rID == 1);
    REQUIRE(gyper->vcf_record.id == "rs0000001");

    REQUIRE(gyper->read_vcf_record() == 1);

    delete gyper;
  }

  SECTION("Same results as above with a compressed VCF file (BCF).")
  {
    std::stringstream base_path;
    base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example_bgz.vcf.gz";
    std::string base_path_str = base_path.str();

    Options CO = Options();
    Gyper* gyper = new Gyper(CO);
    gyper->open_tabix(base_path_str.c_str());

    REQUIRE(gyper->read_tabix_record() == 0);
    std::cout << gyper->tabix_line << std::endl;
    REQUIRE(gyper->read_tabix_record() == 0);
    std::cout << gyper->tabix_line << std::endl;
    REQUIRE(gyper->read_tabix_record() == 0);
    std::cout << gyper->tabix_line << std::endl;
    REQUIRE(gyper->read_tabix_record() == 0);
    std::cout << gyper->tabix_line << std::endl;
    REQUIRE(gyper->read_tabix_record() == 0);
    std::cout << gyper->tabix_line << std::endl;
    REQUIRE(gyper->read_tabix_record() == 1);

    // REQUIRE( == 0);
  //   REQUIRE(gyper->tabix_record.beginPos == 1);
  //   REQUIRE(gyper->tabix_record.id == "rs6054257");
  //   REQUIRE(gyper->tabix_record.ref == "G");
  //   REQUIRE(gyper->tabix_record.alt == "A");
  //   REQUIRE(gyper->tabix_record.qual == 0.0);
  //   REQUIRE(gyper->tabix_record.filter == ".");
  //   REQUIRE(gyper->tabix_record.info == ".");
  //   REQUIRE(gyper->tabix_record.format == "");
  //   REQUIRE(length(gyper->tabix_record.genotypeInfos) == 0);
  //   seqan::String<seqan::Dna> ref_dna = "G";
  //   REQUIRE(static_cast<seqan::String<seqan::Dna>>(gyper->tabix_record.ref) == ref_dna);

  //   REQUIRE(gyper->read_vcf_record() == 0);
  //   REQUIRE(gyper->tabix_record.rID == 0);
  //   REQUIRE(gyper->tabix_record.id == "ins504320");
  //   seqan::String<seqan::Dna> alt_dna = "GT";
  //   REQUIRE(static_cast<seqan::String<seqan::Dna>>(gyper->tabix_record.alt) == alt_dna);

  //   REQUIRE(gyper->read_vcf_record() == 0);
  //   REQUIRE(gyper->tabix_record.rID == 0);
  //   REQUIRE(gyper->tabix_record.id == "del504320");

  //   REQUIRE(gyper->read_vcf_record() == 0);
  //   REQUIRE(gyper->tabix_record.rID == 0);
  //   REQUIRE(gyper->tabix_record.id == "msat14123");

  //   REQUIRE(gyper->read_vcf_record() == 0);
  //   REQUIRE(gyper->tabix_record.rID == 1);
  //   REQUIRE(gyper->tabix_record.id == "rs6054257");

  //   REQUIRE(gyper->read_vcf_record() == 1);

    delete gyper;
  }
}

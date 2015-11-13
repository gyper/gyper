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

// #include <rocksdb/db.h>
// #include <rocksdb/slice.h>
// #include <rocksdb/options.h>


TEST_CASE("BASIC")
{
  SECTION("Pre conditions")
  {
    SECTION("Catch works")
    {
      REQUIRE(true); // yay
    }
  }

  SECTION("Post conditions")
  {
    SECTION("Empty constructor")
    {
      gyper::Constructor* constructor = new gyper::Constructor();
      REQUIRE(constructor != nullptr);
      delete constructor;
    }

    SECTION("Options can be created")
    {
      gyper::Options* CO = new gyper::Options();
      REQUIRE(CO != nullptr);
      delete CO;
    }
    

    SECTION("Specific options")
    {
      gyper::Options* CO = new gyper::Options();
      REQUIRE(CO != nullptr);
      gyper::Constructor* constructor = new gyper::Constructor(*CO);
      REQUIRE(constructor != nullptr);
      delete constructor;
      delete CO;
    }

    SECTION("Setting genomic region")
    {
      gyper::Constructor* constructor = new gyper::Constructor();
      REQUIRE(constructor->genomic_region.beginPos == static_cast<int>(seqan::VcfRecord::INVALID_POS));

      // Read a reference genome, this is tested later
      std::stringstream reference_path;
      reference_path << gyper_SOURCE_DIRECTORY << "/test/reference/reference.fa";
      std::string reference_path_str = reference_path.str();
      constructor->read_reference_genome(reference_path_str.c_str());

      SECTION("chrX:A-B format")
      {
        constructor->set_genomic_region("chr1:2-18");
        REQUIRE(constructor->genomic_region.seqName == "chr1");
        REQUIRE(constructor->genomic_region.rID == 0);
        REQUIRE(constructor->genomic_region.beginPos == 1);
        REQUIRE(constructor->genomic_region.endPos == 18);
      }

      SECTION("chrX:A-B format with a large B")
      {
        constructor->set_genomic_region("chr1:2-999");
        REQUIRE(constructor->genomic_region.seqName == "chr1");
        REQUIRE(constructor->genomic_region.rID == 0);
        REQUIRE(constructor->genomic_region.beginPos == 1);
        REQUIRE(constructor->genomic_region.endPos == 22);
      }

      SECTION("chrX:A format")
      {
        constructor->set_genomic_region("chr1:1");
        REQUIRE(constructor->genomic_region.seqName == "chr1");
        REQUIRE(constructor->genomic_region.rID == 0);
        REQUIRE(constructor->genomic_region.beginPos == 0);
        REQUIRE(constructor->genomic_region.endPos == 22);
      }

      SECTION("chrX format")
      {
        constructor->set_genomic_region("chr2");
        REQUIRE(constructor->genomic_region.seqName == "chr2");
        REQUIRE(constructor->genomic_region.rID == 1);
        REQUIRE(constructor->genomic_region.beginPos == 0);
        REQUIRE(constructor->genomic_region.endPos == 22);
      }
    }
  }
}


TEST_CASE("FASTA I/O")
{
  // Get the base path to test reference genomes
  std::stringstream reference_path;
  reference_path << gyper_SOURCE_DIRECTORY << "/test/reference/reference.fa";
  std::string reference_path_str = reference_path.str();
  gyper::Options* CO = new gyper::Options();
  gyper::Constructor* constructor = new gyper::Constructor(*CO);
  delete CO;

  SECTION("The 'read_reference_genome' method.")
  {
    SECTION("Pre conditions")
    {
      SECTION("Gyper constructor should not pointing to something other than NULL.")
      {
        REQUIRE(constructor != nullptr);
      }

      SECTION("The path to the reference FASTA can be accessed.")
      {
        std::fstream filestr;
        filestr.open(reference_path_str);
        REQUIRE(filestr.is_open());
      }
    }

    SECTION("Post conditions")
    {
      SECTION("All sequences of the FASTA are read.")
      {
        REQUIRE(seqan::numSeqs(constructor->fasta_index) == 0);
        constructor->read_reference_genome(reference_path_str.c_str());
        REQUIRE(seqan::numSeqs(constructor->fasta_index) == 2);
      }
      
      SECTION("We can use the index to query a infix sequence.")
      {
        constructor->read_reference_genome(reference_path_str.c_str());
        unsigned fasta_index_id;

        SECTION("Both 'chr1' and 'chr2' are in the reference FASTA. An nothing else.")
        {
          REQUIRE(seqan::getIdByName(fasta_index_id, constructor->fasta_index, "chr1"));
          REQUIRE(fasta_index_id == 0);
          REQUIRE(seqan::getIdByName(fasta_index_id, constructor->fasta_index, "chr2"));
          REQUIRE(fasta_index_id == 1);
          REQUIRE(seqan::getIdByName(fasta_index_id, constructor->fasta_index, "chr3") == false);
          REQUIRE(fasta_index_id == 1); // unchanged
          REQUIRE(seqan::getIdByName(fasta_index_id, constructor->fasta_index, "chr1"));
          REQUIRE(fasta_index_id == 0); // changed
          REQUIRE(seqan::getIdByName(fasta_index_id, constructor->fasta_index, "asdfas") == false);
          REQUIRE(seqan::getIdByName(fasta_index_id, constructor->fasta_index, "chr") == false);
          REQUIRE(fasta_index_id == 0); // unchanged
        }

        seqan::getIdByName(fasta_index_id, constructor->fasta_index, "chr2");
        seqan::String<seqan::Dna5> sequence_infix;
        seqan::readRegion(sequence_infix, constructor->fasta_index, fasta_index_id, 0u, 1u);
        REQUIRE(sequence_infix == "N");
        seqan::readRegion(sequence_infix, constructor->fasta_index, fasta_index_id, 0u, 100u);
        REQUIRE(sequence_infix == "NNAGGTTTCCCCAGGTTTCCCC");
        seqan::readRegion(sequence_infix, constructor->fasta_index, fasta_index_id, 1u, 0u);
        REQUIRE(sequence_infix == "");
      }
    }
  }

  SECTION("The 'extract_reference_sequence' method.")
  {
    constructor->read_reference_genome(reference_path_str.c_str());

    SECTION("Pre conditions.")
    {
      REQUIRE(seqan::numSeqs(constructor->fasta_index) == 2);
      REQUIRE(seqan::length(constructor->reference_sequence) == 0);
    }

    SECTION("Post conditions.")
    {
      SECTION("We can read an infix of chromosome 2 using the format chrX:A-B.")
      {
        constructor->extract_reference_sequence("chr2:5-10");
        REQUIRE(seqan::length(constructor->reference_sequence) == 6);
        REQUIRE(constructor->reference_sequence == "GTTTCC");
      }

      SECTION("We can read a suffix of chromosome 2 using the format chrX:A.")
      {
        constructor->extract_reference_sequence("chr2:5");
        REQUIRE(seqan::length(constructor->reference_sequence) == 18);
        REQUIRE(constructor->reference_sequence == "GTTTCCCCAGGTTTCCCC");
      }

      SECTION("We can read an entire chromosome 1 using the format chrX.")
      {
        constructor->extract_reference_sequence("chr1");
        REQUIRE(seqan::length(constructor->reference_sequence) == 22);
        REQUIRE(constructor->reference_sequence == "AGGTTTCCCCNNNNNNAGGTTT");
      }

      SECTION("Reading a too large infix will not cause any errors, we simply read as long as we can.")
      {
        constructor->extract_reference_sequence("chr2:5-10000000");
        REQUIRE(seqan::length(constructor->reference_sequence) == 18);
        REQUIRE(constructor->reference_sequence == "GTTTCCCCAGGTTTCCCC");
      }
    }
  }

  delete constructor;
}


TEST_CASE("VCF I/O")
{
  std::stringstream vcf_path;
  vcf_path << gyper_SOURCE_DIRECTORY << "/test/reference/example_bgz.vcf.gz";
  std::string vcf_path_str = vcf_path.str();
  gyper::Options* CO = new gyper::Options();
  gyper::Constructor* constructor = new gyper::Constructor(*CO);
  delete CO;

  SECTION("Opening a Tabix file")
  {
    SECTION("Pre conditions")
    {
      REQUIRE(constructor->vcf_record.rID == static_cast<int>(seqan::VcfRecord::INVALID_REFID));
      REQUIRE(constructor->tabix_file.fp == nullptr);

      std::fstream filestr;
      filestr.open(vcf_path_str);
      REQUIRE(filestr.is_open());
    }

    SECTION("Post conditions")
    {
      SECTION("Entire VCF is read")
      {
        constructor->open_tabix(vcf_path_str.c_str());
        REQUIRE(constructor->tabix_file.fp != nullptr);
        // No record has been read
        REQUIRE(constructor->vcf_record.rID == static_cast<int>(seqan::VcfRecord::INVALID_REFID));
      }

      SECTION("A region is read")
      {
        constructor->open_tabix(vcf_path_str.c_str(), "chr2:2-6");
        REQUIRE(constructor->tabix_file.fp != nullptr);
        REQUIRE(constructor->vcf_record.rID == static_cast<int>(seqan::VcfRecord::INVALID_REFID));
      }
    }
  }

  SECTION("Reading VCF records")
  {
    SECTION("Pre conditions")
    {
      REQUIRE(constructor->vcf_record.rID == static_cast<int>(seqan::VcfRecord::INVALID_REFID));
      REQUIRE(constructor->tabix_file.fp == nullptr);
    }

    SECTION("Post conditions")
    {
      SECTION("Entire VCF is read")
      {
        constructor->open_tabix(vcf_path_str.c_str());
        REQUIRE(constructor->vcf_record.rID == static_cast<int>(seqan::VcfRecord::INVALID_REFID));

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 1);
        REQUIRE(constructor->vcf_record.id == "rs6054257");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 3);
        REQUIRE(constructor->vcf_record.id == "ins504320");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 5);
        REQUIRE(constructor->vcf_record.id == "del504320");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 10);
        REQUIRE(constructor->vcf_record.id == "msat14123");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 1);
        REQUIRE(constructor->vcf_record.beginPos == 2);
        REQUIRE(constructor->vcf_record.id == "rs0000001");

        REQUIRE(!constructor->read_tabix_record());
      }

      SECTION("The VCF is partially read")
      {
        constructor->open_tabix(vcf_path_str.c_str(), "chr1:2-6");
        REQUIRE(constructor->vcf_record.rID == static_cast<int>(seqan::VcfRecord::INVALID_REFID));

        REQUIRE(constructor->read_tabix_region());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 1);
        REQUIRE(constructor->vcf_record.id == "rs6054257");

        REQUIRE(constructor->read_tabix_region());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 3);
        REQUIRE(constructor->vcf_record.id == "ins504320");

        REQUIRE(constructor->read_tabix_region());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 5);
        REQUIRE(constructor->vcf_record.id == "del504320");

        REQUIRE(!constructor->read_tabix_region());
      }
    }
  }

  delete constructor;
}


TEST_CASE("Graph creation")
{
  std::stringstream vcf_path;
  vcf_path << gyper_SOURCE_DIRECTORY << "/test/reference/example6_bgz.vcf.gz";
  std::string vcf_path_str = vcf_path.str();

  std::stringstream reference_path;
  reference_path << gyper_SOURCE_DIRECTORY << "/test/reference/reference.fa";
  std::string reference_path_str = reference_path.str();

  gyper::Options* CO = new gyper::Options();
  gyper::Constructor* constructor = new gyper::Constructor(*CO);
  delete CO;

  SECTION("Pre conditions")
  {
    std::fstream ref_file;
    ref_file.open(reference_path_str);
    REQUIRE(ref_file.is_open());

    std::fstream vcf_file;
    vcf_file.open(vcf_path_str);
    REQUIRE(vcf_file.is_open());
  }

  SECTION("The insertion of vertices")
  {
    SECTION("Pre conditions")
    {
      REQUIRE(seqan::empty(constructor->graph));
      REQUIRE(seqan::numVertices(constructor->graph) == 0);
      REQUIRE(seqan::numEdges(constructor->graph) == 0);
      REQUIRE(constructor->vertex_labels.size() == 0);
    }

    SECTION("Post conditions")
    {
      SECTION("Without a previous vertex")
      {
        gyper::TVertex new_vertex = constructor->insert_reference_vertex(0, "");

        REQUIRE(!seqan::empty(constructor->graph));
        REQUIRE(seqan::numVertices(constructor->graph) == 1);
        REQUIRE(seqan::numEdges(constructor->graph) == 0);
        REQUIRE(seqan::degree(constructor->graph, new_vertex) == 0);

        REQUIRE(constructor->vertex_labels.size() == 1);
        REQUIRE(constructor->vertex_labels[0].order == 0);
        REQUIRE(constructor->vertex_labels[0].dna == "");
      }

      SECTION("Using a previous vertex")
      {
        gyper::TVertex first_vertex = constructor->insert_reference_vertex(0, "");
        gyper::TVertex second_vertex = constructor->insert_reference_vertex(1, "", first_vertex);

        REQUIRE(!seqan::empty(constructor->graph));
        REQUIRE(seqan::numVertices(constructor->graph) == 2);
        REQUIRE(seqan::numEdges(constructor->graph) == 1);

        REQUIRE(seqan::degree(constructor->graph, first_vertex) == 1);
        REQUIRE(seqan::outDegree(constructor->graph, first_vertex) == 1);
        REQUIRE(seqan::inDegree(constructor->graph, first_vertex) == 0);

        REQUIRE(seqan::degree(constructor->graph, second_vertex) == 1);
        REQUIRE(seqan::outDegree(constructor->graph, second_vertex) == 0);
        REQUIRE(seqan::inDegree(constructor->graph, second_vertex) == 1);

        REQUIRE(constructor->vertex_labels.size() == 2);
        REQUIRE(constructor->vertex_labels[0].order == 0);
        REQUIRE(constructor->vertex_labels[0].dna == "");
        REQUIRE(constructor->vertex_labels[1].order == 1);
        REQUIRE(constructor->vertex_labels[1].dna == "");
      }

      SECTION("With some DNA values")
      {
        gyper::TVertex first_vertex = constructor->insert_reference_vertex(0, "A");
        gyper::TVertex second_vertex = constructor->insert_reference_vertex(1, "CT", first_vertex);

        REQUIRE(!seqan::empty(constructor->graph));
        REQUIRE(seqan::numVertices(constructor->graph) == 2);
        REQUIRE(seqan::numEdges(constructor->graph) == 1);

        REQUIRE(seqan::degree(constructor->graph, first_vertex) == 1);
        REQUIRE(seqan::outDegree(constructor->graph, first_vertex) == 1);
        REQUIRE(seqan::inDegree(constructor->graph, first_vertex) == 0);

        REQUIRE(seqan::degree(constructor->graph, second_vertex) == 1);
        REQUIRE(seqan::outDegree(constructor->graph, second_vertex) == 0);
        REQUIRE(seqan::inDegree(constructor->graph, second_vertex) == 1);

        REQUIRE(constructor->vertex_labels.size() == 2);
        REQUIRE(constructor->vertex_labels[0].order == 0);
        REQUIRE(constructor->vertex_labels[0].dna == "A");
        REQUIRE(constructor->vertex_labels[1].order == 1);
        REQUIRE(constructor->vertex_labels[1].dna == "CT");
      }
    }
  }

  SECTION("The 'add_reference_sequence_preceding_a_point' method")
  {
    SECTION("Post conditions")
    {
      constructor->read_reference_genome(reference_path_str.c_str());

      SECTION("Reading a single DNA base from reference")
      {
        constructor->extract_reference_sequence("chr1:1-3");
        const unsigned point = 1; // This is index is 0-based
        constructor->add_first_reference_sequence();
        constructor->add_reference_sequence_preceding_a_point(point);
        REQUIRE(constructor->head == 1);
        REQUIRE(constructor->vertex_labels.size() == 2);
        REQUIRE(constructor->vertex_labels[1].dna == "A");
        REQUIRE(constructor->vertex_labels[1].order == 0);
      }

      SECTION("Reading multiple DNA bases from reference (8 in this case)")
      {
        constructor->extract_reference_sequence("chr1");
        const unsigned point = 8;
        constructor->add_first_reference_sequence();
        constructor->add_reference_sequence_preceding_a_point(point);
        REQUIRE(constructor->head == 1);
        REQUIRE(constructor->vertex_labels.size() == 2);
        REQUIRE(constructor->vertex_labels[1].dna == "AGGTTTCC");
        REQUIRE(constructor->vertex_labels[1].order == 0);
      }

      SECTION("Reading multiple DNA bases with some N's")
      {
        constructor->extract_reference_sequence("chr1");
        const unsigned point = 22; // The entire chr1
        constructor->add_first_reference_sequence();
        constructor->add_reference_sequence_preceding_a_point(point);
        REQUIRE(constructor->head == 3);
        REQUIRE(constructor->vertex_labels.size() == 4);
        REQUIRE(constructor->vertex_labels[0].dna == "");
        REQUIRE(constructor->vertex_labels[0].order == 0);
        REQUIRE(constructor->vertex_labels[1].dna == "AGGTTTCCCC");
        REQUIRE(constructor->vertex_labels[1].order == 0);
        REQUIRE(constructor->vertex_labels[2].dna == "");
        REQUIRE(constructor->vertex_labels[2].order == 16);
        REQUIRE(constructor->vertex_labels[3].dna == "AGGTTT");
        REQUIRE(constructor->vertex_labels[3].order == 16);
      }

      SECTION("Reading a chr that starts with N's")
      {
        constructor->extract_reference_sequence("chr2");
        const unsigned point = 22; // The entire chr2
        constructor->add_first_reference_sequence();
        constructor->add_reference_sequence_preceding_a_point(point);
        REQUIRE(constructor->head == 1);
        REQUIRE(constructor->vertex_labels.size() == 2);
        REQUIRE(constructor->vertex_labels[0].dna == "");
        REQUIRE(constructor->vertex_labels[0].order == 2);
        REQUIRE(constructor->vertex_labels[1].dna == "AGGTTTCCCCAGGTTTCCCC");
        REQUIRE(constructor->vertex_labels[1].order == 2);
      }
    }
  }


  SECTION("The 'add_sequence_preceding_a_vcf_record' method")
  {
    SECTION("Post conditions")
    {
      constructor->read_reference_genome(reference_path_str.c_str());

      SECTION("First VCF record")
      {
        constructor->set_genomic_region("chr1:1-2");
        constructor->open_tabix(vcf_path_str.c_str());
        REQUIRE(constructor->read_first_tabix_region());
        constructor->extract_reference_sequence();
        constructor->add_first_reference_sequence();
        REQUIRE(constructor->head == 0);

        constructor->add_sequence_preceding_a_vcf_record();
        REQUIRE(constructor->head == 1);
        REQUIRE(constructor->vertex_labels.size() == 2);
        REQUIRE(constructor->vertex_labels[0].dna == "");
        REQUIRE(constructor->vertex_labels[0].order == 0);
        REQUIRE(constructor->vertex_labels[1].dna == "A");
        REQUIRE(constructor->vertex_labels[1].order == 0);

        REQUIRE(!constructor->read_tabix_region()); // No more records
      }

      SECTION("Some VCF record")
      {
        constructor->set_genomic_region("chr1:7");
        constructor->open_tabix(vcf_path_str.c_str());
        REQUIRE(constructor->read_first_tabix_region());
        constructor->extract_reference_sequence();
        constructor->add_first_reference_sequence();
        REQUIRE(constructor->head == 0);
        constructor->add_sequence_preceding_a_vcf_record();

        REQUIRE(constructor->head == 3);
        REQUIRE(constructor->vertex_labels.size() == 4);
        REQUIRE(constructor->vertex_labels[0].dna == "");
        REQUIRE(constructor->vertex_labels[0].order == 0);
        REQUIRE(constructor->vertex_labels[1].dna == "CCCC");
        REQUIRE(constructor->vertex_labels[1].order == 0);
        REQUIRE(constructor->vertex_labels[2].dna == "");
        REQUIRE(constructor->vertex_labels[2].order == 10);
        REQUIRE(constructor->vertex_labels[3].dna == "A");
        REQUIRE(constructor->vertex_labels[3].order == 10);

        REQUIRE(!constructor->read_tabix_region()); // No more records
      }

      SECTION("Some other VCF record")
      {
        constructor->set_genomic_region("chr1:8");
        constructor->open_tabix(vcf_path_str.c_str());
        REQUIRE(constructor->read_first_tabix_region());
        constructor->extract_reference_sequence();
        constructor->add_first_reference_sequence();
        REQUIRE(constructor->head == 0);
        constructor->add_sequence_preceding_a_vcf_record();

        REQUIRE(constructor->head == 3);
        REQUIRE(constructor->vertex_labels.size() == 4);
        REQUIRE(constructor->vertex_labels[0].dna == "");
        REQUIRE(constructor->vertex_labels[0].order == 0);
        REQUIRE(constructor->vertex_labels[1].dna == "CCC");
        REQUIRE(constructor->vertex_labels[1].order == 0);
        REQUIRE(constructor->vertex_labels[2].dna == "");
        REQUIRE(constructor->vertex_labels[2].order == 9);
        REQUIRE(constructor->vertex_labels[3].dna == "A");
        REQUIRE(constructor->vertex_labels[3].order == 9);

        REQUIRE(!constructor->read_tabix_region()); // No more records
      }
    } 
  }

  SECTION("Add a VCF record to graph")
  {
    SECTION("Post conditions")
    {
      constructor->read_reference_genome(reference_path_str.c_str());

      SECTION("VCF record with a SNP")
      {
        constructor->set_genomic_region("chr1:1-3");
        constructor->open_tabix(vcf_path_str.c_str());
        REQUIRE(constructor->read_first_tabix_region());
        constructor->extract_reference_sequence();
        constructor->add_first_reference_sequence();
        constructor->add_sequence_preceding_a_vcf_record();

        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 1);
        REQUIRE(constructor->vcf_record.id == "rs6054257");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "A");

        REQUIRE(constructor->vertex_labels.size() == 2);
        constructor->add_vcf_record_to_graph();
        REQUIRE(constructor->vertex_labels.size() == 4);
        REQUIRE(constructor->vertex_labels[2].dna == "G");
        REQUIRE(constructor->vertex_labels[2].order == 1);
        REQUIRE(constructor->vertex_labels[3].dna == "A");
        REQUIRE(constructor->vertex_labels[3].order == 1);

        REQUIRE(!constructor->read_tabix_region()); // No more records
      }

      SECTION("VCF record with two alternative paths")
      {
        constructor->set_genomic_region("chr1:17");
        constructor->open_tabix(vcf_path_str.c_str());
        REQUIRE(constructor->read_first_tabix_region());
        constructor->extract_reference_sequence();
        constructor->add_first_reference_sequence();
        constructor->add_sequence_preceding_a_vcf_record();

        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 17);
        REQUIRE(constructor->vcf_record.id == "msat14123");
        REQUIRE(constructor->vcf_record.ref == "GGT");
        REQUIRE(constructor->vcf_record.alt == "G,GGTGTGT");

        REQUIRE(constructor->vertex_labels.size() == 2);
        constructor->add_vcf_record_to_graph();
        REQUIRE(constructor->vertex_labels.size() == 5);
        REQUIRE(constructor->vertex_labels[2].dna == "GGT");
        REQUIRE(constructor->vertex_labels[2].order == 17);
        REQUIRE(constructor->vertex_labels[3].dna == "G");
        REQUIRE(constructor->vertex_labels[3].order == 17);
        REQUIRE(constructor->vertex_labels[4].dna == "GGTGTGT");
        REQUIRE(constructor->vertex_labels[4].order == 17);

        REQUIRE(!constructor->read_tabix_region()); // No more records
      }
    }
  }

  delete constructor;
}

TEST_CASE("THE ONLY THING YOU NEED")
{
  std::stringstream vcf_path;
  vcf_path << gyper_SOURCE_DIRECTORY << "/test/reference/example6_bgz.vcf.gz";
  std::string vcf_path_str = vcf_path.str();

  std::stringstream reference_path;
  reference_path << gyper_SOURCE_DIRECTORY << "/test/reference/reference.fa";
  std::string reference_path_str = reference_path.str();

  gyper::Options* CO = new gyper::Options();
  gyper::Constructor* constructor = new gyper::Constructor(*CO);
  delete CO;

  std::string expected = "";

  SECTION("chr1")
  {
    constructor->construct_graph(reference_path_str.c_str(), vcf_path_str.c_str(), "chr1");
    expected = std::string("Head is at: 16\nVertex 0:  @ 0\nVertex 1: A @ 0\nVertex 2: G @ 1\nVertex 3: A @ 1\nVertex 4: G @ 2\nVertex 5: G @ 3\nVertex 6: GT @ 3\nVertex 7: T @ 4\nVertex 8: TC @ 5\nVertex 9: T @ 5\nVertex 10: CCC @ 7\nVertex 11:  @ 16\nVertex 12: A @ 16\nVertex 13: GGT @ 17\nVertex 14: G @ 17\nVertex 15: GGTGTGT @ 17\nVertex 16: TT @ 20\n");
  }

  SECTION("chr2")
  {
    constructor->construct_graph(reference_path_str.c_str(), vcf_path_str.c_str(), "chr2");

    expected = std::string("Head is at: 3\nVertex 0:  @ 2\nVertex 1: A @ 2\nVertex 2: T @ 2\nVertex 3: GGTTTCCCCAGGTTTCCCC @ 3\n");
// Vertex 1: A @ 0
// Vertex 2: G @ 1
// Vertex 3: A @ 1
// Vertex 4: G @ 2
// Vertex 5: G @ 3
// Vertex 6: GT @ 3
// Vertex 7: T @ 4
// Vertex 8: TC @ 5
// Vertex 9: T @ 5
// Vertex 10: CCC @ 7
// Vertex 11:  @ 16
// Vertex 12: A @ 16
// Vertex 13: GGT @ 17
// Vertex 14: G @ 17
// Vertex 15: GGTGTGT @ 17
// Vertex 16: TT @ 20

  }

  std::ostringstream ss;
  ss << *constructor;

  REQUIRE(ss.str() == expected);
  // std::cout << "constructor = " << *constructor << std::endl;
}

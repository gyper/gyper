#include <stdio.h>
#include <climits>

#include "../src/partial_order_graph.hpp"
#include "../src/gyper_options.hpp"
// #include "../src/graph.hpp"
#include <catch.hpp>
#include "constants.hpp"

// #include "rocksdb/db.h"
// #include "rocksdb/slice.h"
// #include "rocksdb/options.h"

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
		Gyper* gyper = new gyper::Gyper(CO);
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
		REQUIRE(gyper->vcf_record.ref == "G");
		REQUIRE(gyper->vcf_record.alt == "GT");
		seqan::String<seqan::Dna> alt_dna = "GT";
		REQUIRE(static_cast<seqan::String<seqan::Dna>>(gyper->vcf_record.alt) == alt_dna);

		REQUIRE(gyper->read_vcf_record() == 0);
		REQUIRE(gyper->vcf_record.rID == 0);
		REQUIRE(gyper->vcf_record.id == "del504320");
		REQUIRE(gyper->vcf_record.ref == "TC");
		REQUIRE(gyper->vcf_record.alt == "T");

		REQUIRE(gyper->read_vcf_record() == 0);
		REQUIRE(gyper->vcf_record.rID == 0);
		REQUIRE(gyper->vcf_record.id == "msat14123");
		REQUIRE(gyper->vcf_record.ref == "GGT");
		REQUIRE(gyper->vcf_record.alt == "G,GGTGTGT");

		REQUIRE(gyper->read_vcf_record() == 0);
		REQUIRE(gyper->vcf_record.rID == 1);
		REQUIRE(gyper->vcf_record.id == "rs0000001");
		REQUIRE(gyper->vcf_record.ref == "A");
		REQUIRE(gyper->vcf_record.alt == "T");

		REQUIRE(gyper->read_vcf_record() == 1);

		delete gyper;
	}

	SECTION("BCF (binary vcf)")
	{
		SECTION("Same results as above with a compressed VCF file (BCF) (example_bgz.vcf.gz).")
		{
			std::stringstream base_path;
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example_bgz.vcf.gz";
			std::string base_path_str = base_path.str();

			Options CO = Options();
			Gyper* gyper = new gyper::Gyper(CO);
			gyper->open_tabix(base_path_str.c_str());

			SECTION("The getNextFunction")
			{
				// REQUIRE(gyper->read_tabix_record() == 0);
				seqan::String<char> line;
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	2	rs6054257	G	A	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	4	ins504320	G	GT	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	6	del504320	TC	T	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	11	msat14123	GGT	G,GGTGTGT	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr2	3	rs0000001	A	T	0	.	.");
				REQUIRE(!seqan::readRawRecord(line, gyper->tabix_file));
			}

			SECTION("The readRecord method")
			{
				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.beginPos == 1);
				REQUIRE(gyper->vcf_record.id == "rs6054257");
				REQUIRE(gyper->vcf_record.ref == "G");
				REQUIRE(gyper->vcf_record.alt == "A");
				REQUIRE(gyper->vcf_record.qual == 0.0);
				REQUIRE(gyper->vcf_record.filter == ".");
				REQUIRE(gyper->vcf_record.info == ".");
				REQUIRE(gyper->vcf_record.format == "");
				REQUIRE(seqan::length(gyper->vcf_record.genotypeInfos) == 0);
				seqan::String<seqan::Dna> ref_dna = "G";
				REQUIRE(static_cast<seqan::String<seqan::Dna>>(gyper->vcf_record.ref) == ref_dna);

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "ins504320");
				REQUIRE(gyper->vcf_record.ref == "G");
				REQUIRE(gyper->vcf_record.alt == "GT");
				seqan::String<seqan::Dna> alt_dna = "GT";
				REQUIRE(static_cast<seqan::String<seqan::Dna>>(gyper->vcf_record.alt) == alt_dna);

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "del504320");
				REQUIRE(gyper->vcf_record.ref == "TC");
				REQUIRE(gyper->vcf_record.alt == "T");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "msat14123");
				REQUIRE(gyper->vcf_record.ref == "GGT");
				REQUIRE(gyper->vcf_record.alt == "G,GGTGTGT");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 1);
				REQUIRE(gyper->vcf_record.id == "rs0000001");
				REQUIRE(gyper->vcf_record.ref == "A");
				REQUIRE(gyper->vcf_record.alt == "T");

				REQUIRE(!gyper->read_tabix_record());
			}

			SECTION("The readRegion method")
			{
				seqan::String<seqan::VcfRecord> records;

				SECTION("Using extractions of all variants in whole chromosomes")
				{
					// Read chr1
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1");
					REQUIRE(seqan::length(records) == 4);

					// Read chr2
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr2");
					REQUIRE(seqan::length(records) == 1);

					// Read both chromosomes
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1");
					readRegion(records, gyper->tabix_file, "chr2");
					REQUIRE(seqan::length(records) == 5);
				}
				
				SECTION("Using extraction of partial chromosomes")
				{
					// Read chr1:2-2. This should only read the variant on position 2 of chr 1.
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1:2-2");
					REQUIRE(seqan::length(records) == 1);

					// Read chr1:2-3. This does not read the variant on position 4.
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1:2-3");
					REQUIRE(seqan::length(records) == 1);

					// Read chr1:2-3
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1:2-4");
					REQUIRE(seqan::length(records) == 2);

					// Read chr1:4-4
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1:4-4");
					REQUIRE(seqan::length(records) == 1);

					// Read chr1:4-6
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1:4-6");
					REQUIRE(seqan::length(records) == 2);

					// Read chr1:0-10000000000. This reads the entire chromosome 1
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr1:1-10000000000");
					REQUIRE(seqan::length(records) == 4);

					// Read chr2:0-0.
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr2:0-0");
					REQUIRE(seqan::length(records) == 0);

					// Read chr3:3-3.
					clear(records);
					REQUIRE(seqan::length(records) == 0);
					readRegion(records, gyper->tabix_file, "chr2:3-3");
					REQUIRE(seqan::length(records) == 1);
				}
			}

			delete gyper;
		}

		SECTION("Another BCF example with 2 variants on chr1 and 3 on chr2. (example2_bgz.vcf.gz)")
		{
			std::stringstream base_path;
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example2_bgz.vcf.gz";
			std::string base_path_str = base_path.str();

			Options CO = Options();
			Gyper* gyper = new gyper::Gyper(CO);
			gyper->open_tabix(base_path_str.c_str());

			SECTION("The getNextFunction")
			{
				seqan::String<char> line;
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	2	rs6054257	G	A	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	4	ins504320	T	GT	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr2	6	del504320	T	TA	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr2	13	del0000001	AG	A	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr2	14	mix14123	G	GTT,T	0	.	.");
				REQUIRE(!seqan::readRawRecord(line, gyper->tabix_file));
			}

			SECTION("The readRecord method")
			{
				REQUIRE(gyper->read_tabix_record());
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

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "ins504320");
				REQUIRE(gyper->vcf_record.ref == "T");
				REQUIRE(gyper->vcf_record.alt == "GT");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 1);
				REQUIRE(gyper->vcf_record.id == "del504320");
				REQUIRE(gyper->vcf_record.ref == "T");
				REQUIRE(gyper->vcf_record.alt == "TA");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 1);
				REQUIRE(gyper->vcf_record.id == "del0000001");
				REQUIRE(gyper->vcf_record.ref == "AG");
				REQUIRE(gyper->vcf_record.alt == "A");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 1);
				REQUIRE(gyper->vcf_record.id == "mix14123");
				REQUIRE(gyper->vcf_record.ref == "G");
				REQUIRE(gyper->vcf_record.alt == "GTT,T");

				REQUIRE(!gyper->read_tabix_record());
			}

			delete gyper;
		}

		SECTION("Another BCF example with 2 variants on chr1. (example3.vcf.gz)")
		{
			std::stringstream base_path;
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example3.vcf.gz";
			std::string base_path_str = base_path.str();

			Options CO = Options();
			Gyper* gyper = new gyper::Gyper(CO);
			gyper->open_tabix(base_path_str.c_str());

			SECTION("The getNextFunction")
			{
				seqan::String<char> line;
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	2	rs6054257	G	A	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	4	ins504320	T	GT	0	.	.");
				REQUIRE(!seqan::readRawRecord(line, gyper->tabix_file));
			}

			SECTION("The readRecord method")
			{
				REQUIRE(gyper->read_tabix_record());
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

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "ins504320");
				REQUIRE(gyper->vcf_record.ref == "T");
				REQUIRE(gyper->vcf_record.alt == "GT");

				REQUIRE(!gyper->read_tabix_record());
			}

			delete gyper;
		}


		SECTION("Another BCF example with and 3 on chr2. (example4.vcf.gz)")
		{
			std::stringstream base_path;
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example4.vcf.gz";
			std::string base_path_str = base_path.str();

			Options CO = Options();
			Gyper* gyper = new gyper::Gyper(CO);
			gyper->open_tabix(base_path_str.c_str());

			SECTION("The getNextFunction")
			{
				seqan::String<char> line;
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr2	6	del504320	T	TA	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr2	13	del0000001	AG	A	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr2	14	mix14123	G	GTT,T	0	.	.");
				REQUIRE(!seqan::readRawRecord(line, gyper->tabix_file));
			}

			SECTION("The readRecord method")
			{
				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "del504320");
				REQUIRE(gyper->vcf_record.ref == "T");
				REQUIRE(gyper->vcf_record.alt == "TA");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "del0000001");
				REQUIRE(gyper->vcf_record.ref == "AG");
				REQUIRE(gyper->vcf_record.alt == "A");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.id == "mix14123");
				REQUIRE(gyper->vcf_record.ref == "G");
				REQUIRE(gyper->vcf_record.alt == "GTT,T");

				REQUIRE(!gyper->read_tabix_record());
			}

			delete gyper;
		}

		SECTION("Another BCF example with variants on chr1 and chr3. (example5.vcf.gz)")
		{
			std::stringstream base_path;
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example5.vcf.gz";
			std::string base_path_str = base_path.str();

			Options CO = Options();
			Gyper* gyper = new gyper::Gyper(CO);
			gyper->open_tabix(base_path_str.c_str());

			SECTION("The getNextFunction")
			{
				seqan::String<char> line;
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr1	6	del504320	T	TA	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr3	13	del0000001	AG	A	0	.	.");
				REQUIRE(seqan::readRawRecord(line, gyper->tabix_file));
				REQUIRE(line == "chr3	14	mix14123	G	GTT,T	0	.	.");
				REQUIRE(!seqan::readRawRecord(line, gyper->tabix_file));
			}

			SECTION("The readRecord method")
			{
				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 0);
				REQUIRE(gyper->vcf_record.beginPos == 5);
				REQUIRE(gyper->vcf_record.id == "del504320");
				REQUIRE(gyper->vcf_record.ref == "T");
				REQUIRE(gyper->vcf_record.alt == "TA");
				REQUIRE(gyper->vcf_record.qual == 0.0);
				REQUIRE(gyper->vcf_record.filter == ".");
				REQUIRE(gyper->vcf_record.info == ".");
				REQUIRE(gyper->vcf_record.format == "");
				REQUIRE(length(gyper->vcf_record.genotypeInfos) == 0);

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 1);
				REQUIRE(gyper->vcf_record.beginPos == 12);
				REQUIRE(gyper->vcf_record.id == "del0000001");
				REQUIRE(gyper->vcf_record.ref == "AG");
				REQUIRE(gyper->vcf_record.alt == "A");

				REQUIRE(gyper->read_tabix_record());
				REQUIRE(gyper->vcf_record.rID == 1);
				REQUIRE(gyper->vcf_record.beginPos == 13);
				REQUIRE(gyper->vcf_record.id == "mix14123");
				REQUIRE(gyper->vcf_record.ref == "G");
				REQUIRE(gyper->vcf_record.alt == "GTT,T");

				REQUIRE(!gyper->read_tabix_record());
			}

			delete gyper;
		}
	}
}

TEST_CASE("CRAM support")
{
	SECTION("Sample file test")
	{
		std::stringstream base_path;

		SECTION("SAM")
		{
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.sam";
		}

		SECTION("BAM")
		{
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.bam";
		}

		SECTION("CRAM")
		{
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.cram";
		}
		
		std::string base_path_str = base_path.str();
		seqan::HtsSequenceRecord hts_record;
		seqan::HtsFile hts_file(base_path_str.c_str());
		seqan::open(hts_file);

		while (seqan::readRecord(hts_record, hts_file))
		{
			// std::cout << hts_record.qName << " " << hts_record.seq << std::endl;
		}
	}

	SECTION("HTS indexing")
	{
		std::stringstream base_path;

		SECTION("BAM")
		{
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.bam";
		}

		SECTION("CRAM")
		{
			base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.cram";
		}

		std::string base_path_str = base_path.str();
		seqan::HtsSequenceRecord hts_record;
		seqan::HtsFile hts_file(base_path_str.c_str());
		seqan::open(hts_file);
		seqan::loadIndex(hts_file, true);
		seqan::setRegion(hts_file, "chr1:238-239");

		while (seqan::readRegion(hts_record, hts_file))
		{
			// std::cout << hts_record.qName << " " << hts_record.seq << std::endl;
		}
  }
}

// TEST_CASE("rocksdb test")
// {
// 	std::string kDBPath = "/tmp/rocksdb_simple_example";

// 	rocksdb::DB* db;
//   rocksdb::Options options;
//   // Optimize RocksDB. This is the easiest way to get RocksDB to perform well
//   options.IncreaseParallelism();
//   options.OptimizeLevelStyleCompaction();
//   // create the DB if it's not already present
//   options.create_if_missing = true;

//   // open DB
//   rocksdb::Status s = DB::Open(options, kDBPath, &db);
// }

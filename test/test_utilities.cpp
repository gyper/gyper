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
#include <gyper/utilities.hpp>


TEST_CASE("To uint_64")
{
  SECTION("Empty string")
  {
    seqan::String<seqan::Dna> dna = "";
    REQUIRE(to_uint64(dna) == 0);
    REQUIRE(to_dna(to_uint64(dna), 0) == dna);
  }

  SECTION("One A base")
  {
    seqan::String<seqan::Dna> dna = "A";
    REQUIRE(to_uint64(dna) == 0);
    REQUIRE(to_dna(to_uint64(dna), 1) == dna);
  }

  SECTION("One C base")
  {
    seqan::String<seqan::Dna> dna = "C";
    REQUIRE(to_uint64(dna) == 1);
    REQUIRE(to_dna(to_uint64(dna), 1) == dna);
  }

  SECTION("One C base")
  {
    seqan::String<seqan::Dna> dna = "G";
    REQUIRE(to_uint64(dna) == 2);
    REQUIRE(to_dna(to_uint64(dna), 1) == dna);
  }

  SECTION("One T base")
  {
    seqan::String<seqan::Dna> dna = "T";
    REQUIRE(to_uint64(dna) == 3);
    REQUIRE(to_dna(to_uint64(dna), 1) == dna);
  }

  SECTION("AA")
  {
    seqan::String<seqan::Dna> dna = "AA";
    REQUIRE(to_uint64(dna) == 0);
    REQUIRE(to_dna(to_uint64(dna), 2) == dna);
  }

  SECTION("TT")
  {
    seqan::String<seqan::Dna> dna = "TT";
    REQUIRE(to_uint64(dna) == 15);   // 0000..001111
    REQUIRE(to_dna(to_uint64(dna), 2) == dna);
  }

  SECTION("CTT")
  {
    seqan::String<seqan::Dna> dna = "CTT";
    REQUIRE(to_uint64(dna) == 31);   // 0000..011111
    REQUIRE(to_dna(to_uint64(dna), 3) == dna);
  }

  SECTION("GTT")
  {
    seqan::String<seqan::Dna> dna = "GTT";
    REQUIRE(to_uint64(dna) == 15+32); // 0000..0101111
    REQUIRE(to_dna(to_uint64(dna), 3) == dna);
  }

  SECTION("TTT")
  {
    seqan::String<seqan::Dna> dna = "TTT";
    REQUIRE(to_uint64(dna) == 0x000000000000003F); // 0000..0111111
    REQUIRE(to_dna(to_uint64(dna), 3) == dna);
  }

  SECTION("TTTT")
  {
    seqan::String<seqan::Dna> dna = "TTTT";
    REQUIRE(to_uint64(dna) == 0x00000000000000FF); // 000..01111111
    REQUIRE(to_dna(to_uint64(dna), 4) == dna);
  }

  SECTION("TTTTTTTTTTTTTTTT")
  {
    seqan::String<seqan::Dna> dna = "TTTTTTTTTTTTTTTT";
    REQUIRE(to_uint64(dna) == 0x00000000FFFFFFFF); // 00..0011..11
    REQUIRE(to_dna(to_uint64(dna), 16) == dna);
  }

  SECTION("TTTTTTTTTTTTTTTTT")
  {
    seqan::String<seqan::Dna> dna = "TTTTTTTTTTTTTTTTT";
    REQUIRE(to_uint64(dna) == 0x00000003FFFFFFFF);
    REQUIRE(to_dna(to_uint64(dna), 17) == dna);
  }

  SECTION("32nd DNA base matters")
  {
    seqan::String<seqan::Dna> mer32 = "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    seqan::String<seqan::Dna> mer31 = "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    REQUIRE(to_uint64(mer32) > to_uint64(mer31));
    REQUIRE(to_uint64(mer32) == (static_cast<uint64_t>(1) << 62));
    REQUIRE(to_uint64(mer31) == (static_cast<uint64_t>(1) << 60));
    mer31 = "GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    REQUIRE(to_uint64(mer31) == (static_cast<uint64_t>(1) << 61));

    REQUIRE(to_dna(to_uint64(mer32), 32) == mer32);
    REQUIRE(to_dna(to_uint64(mer31), 31) == mer31);
  }
}

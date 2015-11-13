#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

// functions to split a string by a specific delimiter
#include <stdio.h>
#include <iostream>
#include <stdint.h> // uint64_t

#include <seqan/basic.h>
#include <seqan/sequence.h>


static const uint64_t K_SIZE = 32;                     /** \brief Default size of k is 32 (which fills up all 64 bits). */
static const uint64_t A_VALUE = 0x0000000000000000ULL; /** \brief 'A' is represented by '00'. */
static const uint64_t C_VALUE = 0x0000000000000001ULL; /** \brief 'C' is represented by '01'. */
static const uint64_t G_VALUE = 0x0000000000000002ULL; /** \brief 'G' is represented by '10'. */
static const uint64_t T_VALUE = 0x0000000000000003ULL; /** \brief 'T' is represented by '11'. */

/**
 * @brief Converts a single DNA base to unsigned 64 bit integer.
 * 
 * @param d The DNA base to be converted.
 * @return The new unsigned 64 bit integer.
 */
uint64_t inline
to_uint64(seqan::Dna d)
{
  return seqan::ordValue(d);
}

/**
 * @brief Converts a string/list of DNA bases to a unsigned 64 bit integer.
 * 
 * @param s String of DNA bases.
 * @return The new unsigned 64 bit integer.
 */
uint64_t inline
to_uint64(seqan::String<seqan::Dna> s)
{
  SEQAN_ASSERT_MSG(seqan::length(s) <= 32, "Maximum allowed string length is 32, your string had length %u", seqan::length(s));
  uint64_t d = 0;

  for (unsigned i = 0 ; i < seqan::length(s) ; ++i)
  {
    d *= 4;
    d += seqan::ordValue(s[i]);
  }

  return d;
}

seqan::String<seqan::Dna> inline
to_dna(uint64_t const & d, uint64_t k)
{
  seqan::String<seqan::Dna> new_dna_string = "";

  while (k > 0)
  {
    --k;

    switch ((d & (0x0000000000000003ULL << 2*k)) >> 2*k)
    {
      case A_VALUE:
        seqan::append(new_dna_string, seqan::Dna('A'));
        break;

      case C_VALUE:
        seqan::append(new_dna_string, seqan::Dna('C'));
        break;

      case G_VALUE:
        seqan::append(new_dna_string, seqan::Dna('G'));
        break;

      case T_VALUE:
        seqan::append(new_dna_string, seqan::Dna('T'));
        break;
    }
  }

  return new_dna_string;
}

seqan::String<seqan::Dna> inline
to_dna(uint64_t const & d)
{
  return to_dna(d, K_SIZE);
}

#endif // __UTILITIES_HPP__

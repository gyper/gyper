#include "fasta_region.hpp"

#include <assert.h>
#include <seqan/seq_io.h>

namespace gyper
{

FastaRegion::FastaRegion(unsigned index_id, unsigned begin_position, unsigned end_position)
{
  assert(begin_position <= end_position);
  this->index_id = index_id;
  this->begin_position = begin_position;
  this->end_position = end_position;
}

} // namespace gyper

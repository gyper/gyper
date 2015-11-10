#ifndef __FASTA_REGION_HPP__
#define __FASTA_REGION_HPP__

namespace gyper
{

class FastaRegion
{
 public:
  unsigned index_id;
  unsigned begin_position;
  unsigned end_position;

  FastaRegion(unsigned index_id, unsigned begin_position, unsigned end_position);
};

} // namespace gyper
#endif

#ifndef __GYPER_LABEL_HPP__
#define __GYPER_LABEL_HPP__

#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered/unordered_map.hpp>

#include <gyper/constructor.hpp>

namespace gyper
{

class KmerLabel
{
 public:
  uint32_t order;
  // TVertex start_vertex;
  uint32_t start_index;
  // TVertex end_vertex;
  uint32_t end_index;
  boost::dynamic_bitset<> variant_bits;

  KmerLabel(uint32_t const & order, uint32_t const & start_index, uint32_t const & end_index)
    : order(order), start_index(start_index), end_index(end_index), variant_bits() {}
  
};

// #ifdef GYPER_TESTING
// typedef std::unordered_map< uint8_t, std::forward_list<KmerLabel> > TKmerMap;
// #else // GYPER_TESTING
// typedef std::unordered_map< uint64_t, std::forward_list<KmerLabel> > TKmerMap;
typedef std::unordered_map< uint64_t, std::forward_list<int> > TKmerMap; // testing
// #endif // GYPER_TESTING

} // namespace gyper

#endif // __GYPER_LABEL_HPP__

#ifndef __GYPER_INDEXER_HPP__
#define __GYPER_INDEXER_HPP__

#include <cstdint>
#include <unordered_map>
#include <list>

#include <gyper/constructor.hpp>
#include <gyper/kmer_label.hpp>

namespace gyper
{

class Indexer
{
 public:
  static const uint64_t k = 32;

  TKmerMap kmers;               /** \brief Hashmap with a k-mer as a key and a list of genotyping information as value. */
  Constructor* c;               /** \brief Reference to a constructor which contains a graph. */
  
  unsigned current_mer;
  // #ifdef GYPER_TESTING
  // uint64_t k = 4;
  // #else
  
  // #endif

  std::list<std::deque<uint64_t> > mers; /** \brief A list of short sequences. */
  
  unsigned number_of_upcoming_nodes;


  /***************
   * FIRST STEPS *
   ***************/
  bool is_N_vertex(TVertex const & v);
  void insert_base_to_list(std::list<std::deque<uint64_t> > & list, seqan::Dna const & base);
  void insert_to_map(std::deque<uint64_t> const & q);
  void index_variant(TVertex & vertex);
  void index_reference(TVertex const & vertex);
  void insert_list_into_the_mers_list(std::list<std::deque<uint64_t> > & list);
  void initial_indexing();
  void index_graph();

  Indexer(Constructor* c)
    : kmers(), c(c), current_mer(1), mers(0), number_of_upcoming_nodes(1)
  {

  }

};

} // namespace gyper

#endif // __GYPER_INDEXER_HPP__

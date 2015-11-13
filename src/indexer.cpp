#include <gyper/indexer.hpp>
#include <gyper/utilities.hpp>

// std::deque<uint64_t>
// clear_all(std::deque<uint64_t> q)
// {
//   q.clear();
//   return q;
// }

void
print_all(std::list<std::deque<uint64_t> > q)
{
  std::cout << "List of size " << std::setw(2) << q.size() << " had queues: ";
  unsigned i = 1;

  for (auto q_it = q.begin(); q_it != q.end(); ++q_it, ++i)
  {
    for (auto dq_it = q_it->begin(); dq_it != q_it->end(); ++dq_it)
    {
      std::cout << to_dna(*dq_it, i);
    }

    std::cout << " | ";
  }

  std::cout << std::endl;
}

// std::deque<uint64_t>
// add_(std::deque<uint64_t> const & d)
// {
//   return d;
// }

namespace gyper
{

void
Indexer::insert_list_into_the_mers_list(std::list<std::deque<uint64_t> > & list)
{
  SEQAN_ASSERT_MSG(list.size() == mers.size(), "Both lists should be of the same size.");
  auto mer_it = mers.begin();
  auto list_it = list.begin();

  for (; list_it != mers.end(); ++list_it, ++mer_it)
  {
    for (auto list_kmer_it = list_it->begin(); list_kmer_it != list_it->end(); ++list_kmer_it)
    {
      mer_it->push_back(*list_kmer_it);
    }
  }
}


void
Indexer::index_variant(TVertex & v)
{
  std::list<std::deque<uint64_t> > clean_list(mers); // copies all mers, we find new kmers using the copy.
  // // mers.shift(-1); // Shift to the right

  // Do stuff

  // // Loops over variants
  while (number_of_upcoming_nodes > 0)
  {
    std::cout << "This is a var at " << v << ". " << number_of_upcoming_nodes-1 << " alts remain." << std::endl;
    std::list<std::deque<uint64_t> > new_list(clean_list); // copies all mers, we find new kmers using the copy.



  //   SEQAN_ASSERT(seqan::length(c->vertex_labels[v].dna) != 0); // Empty variant not allowed

  //   std::cout << c->vertex_labels[v].dna << std::endl;

  //   // Loops over DNA bases in the current variant
  //   for (unsigned d = 0; d < seqan::length(c->vertex_labels[v].dna); ++d)
  //   {
  //     // switch()

  //     // Loops over K-mer sizes
  //     // for (unsigned i = 0; i < current_mer; ++i)
  //     // {

  //     // }
  //   }

    --number_of_upcoming_nodes;
    ++v;
  }

  number_of_upcoming_nodes = 1; // Always a reference expected next
}

void
Indexer::insert_to_map(std::deque<uint64_t> const & q)
{
  (void) q;
  std::cout << "Herp derp, inserting to map." << std::endl; 
}

bool inline
Indexer::is_N_vertex(TVertex const & v)
{
  return seqan::length(c->vertex_labels[v].dna) == 0;
}

void
Indexer::insert_base_to_list(std::list<std::deque<uint64_t> > & list, seqan::Dna const & base)
{
  SEQAN_ASSERT_MSG(list.size() < k, "Please remove items from the list before inserting.");
  // Add the base to all existing k-mers
  for (auto mer_it = list.begin(); mer_it != list.end(); ++mer_it)
  {
    for (auto uint64_it = mer_it->begin(); uint64_it != mer_it->end(); ++uint64_it)
    {
      (*uint64_it) <<= 2; // This is exactly the same as multiplying by four, but this is way cooler.
      (*uint64_it) += to_uint64(base);
    }
  }

  // Add a new element with the new DNA base
  std::deque<uint64_t> new_deque(1, to_uint64(base));
  list.push_front(new_deque);

  // std::cout << c->vertex_labels[v].dna << std::endl;
  print_all(list);
}


void
Indexer::index_reference(TVertex const & v)
{
  std::cout << "This is a ref " << v << std::endl;

  if (is_N_vertex(v))
  {
    // Sequences of N's, so we start over with a fresh list. Until next time!
    mers.clear();
  }
  else
  {
    // Loops over any additional DNA bases on the reference
    for (unsigned d = 0; d < seqan::length(c->vertex_labels[v].dna); ++d)
    {
      if (mers.size() >= k)
      {
        insert_to_map(mers.back());
        mers.pop_back();
      }

      insert_base_to_list(mers, c->vertex_labels[v].dna[d]);
    }
  }
}

void
Indexer::index_graph()
{
  TVertex v = 0;

  while (v < seqan::numVertices(c->graph))
  {
    switch (number_of_upcoming_nodes)
    {
      case 1:
        index_reference(v);
        number_of_upcoming_nodes = seqan::outDegree(c->graph, v);
        ++v;
        break;

      default:
        index_variant(v);
    }
  }
}

}


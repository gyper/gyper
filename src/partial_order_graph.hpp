#ifndef __GYPER_HPP__
#define __GYPER_HPP__
#include "graph.hpp"
#include "gyper_options.hpp"

/**
 * @brief Gyper's main class
 * @details The one class to rule them all. Each process should only have one Gyper instance.
 */

class Gyper
{
 private:
  /** @brief Vector of all nodes (vertex) labels. */
  std::vector<VertexLabels> vertex_vector;

  /** @brief Holds a list of labels on each vertex. */
  std::vector<std::string> ids;

  /** @brief Holds a list of labels on each vertex. */
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  /** @brief A set of all nodes with no DNA base. */
  boost::unordered_set<TVertexDescriptor> free_nodes;

  /** @brief The topological order of all nodes (vertexes). */
  String<TVertexDescriptor> order;

  TVertexDescriptor begin_vertex;
  TVertexDescriptor new_begin_vertex;

 public:
  /**
   * @brief Default constructor
   * @details [long description]
   * @return [description]
   * 
   * @post A Gyper object with a empty partial order graph
   */
  Gyper ();

  /**
   * @brief Constructor with specific call options.
   * @details [long description]
   * 
   * @param CO Call options.
   * @return A new Gyper object.
   */
  Gyper (Options & CO);

  /**
   * @brief Creates a single HLA graph.
   * @details Creates an partial order graph which represents a HLA gene specified by the call option.
   */
  void create_HLA_graph();

  /**
   * @brief Indexes the partial order graph
   * @details [long description]
   */
  void index();

  /** @brief graph is a SeqAn graph */
  TGraph graph;

  /** @brief The argument options */
  Options CO;


 private:
  /**
   * @brief Get the number of exons for the selected gene.
   * @details [long description]
   * @return The number of exons this HLA gene has.
   */
  unsigned get_number_of_exons();

  std::string get_HLA_base_path();

  void add_HLA_intron();

  void add_FASTA_region(bool add_bitstrings, int feature_number = 0, bool intron_region = false, bool p3_region = false, bool p5_region = false);
};



#endif

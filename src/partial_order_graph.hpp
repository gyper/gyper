#ifndef __GYPER_HPP__
#define __GYPER_HPP__
// #include "graph.hpp"
#include "gyper_options.hpp"
#include "constants.hpp"
#include "graph_creation.hpp"
#include "fasta_region.hpp"

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>
#include <seqan/bam_io.h>
#include <seqan/stream.h>

#include <iostream>
#include <unordered_map>
#include <boost/unordered/unordered_map.hpp>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>

typedef seqan::Graph<seqan::Directed<void, seqan::WithoutEdgeId> > TGraph;
typedef seqan::VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef seqan::EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

struct KmerLabels
{
  TVertexDescriptor start_vertex;
  TVertexDescriptor end_vertex;
  boost::dynamic_bitset<> id_bits;
};

typedef boost::unordered_map<seqan::String<seqan::Dna>, std::vector<KmerLabels> > TKmerMap;

/**
 * @brief Gyper's main class
 * @details The one class to rule them all. Each process should only have one Gyper instance.
 */

class Gyper
{
 private:
  std::vector<VertexLabel> vertex_labels;

  /** @brief Vector of all nodes (vertex) labels. */
  // std::vector<VertexLabels> vertex_vector;

  /** @brief Holds a list of labels on each vertex. */
  std::vector<std::string> ids;

  /** @brief Holds a list of labels on each vertex. */
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;

  /** @brief A set of all nodes with no DNA base. */
  std::unordered_set<TVertexDescriptor> free_nodes;

  /** @brief The topological order of all nodes (vertexes). */
  seqan::String<TVertexDescriptor> order;

  TVertexDescriptor begin_vertex;

  TVertexDescriptor new_begin_vertex;

  std::map<VertexLabel, TVertexDescriptor> vertex_label_map;


 public:
   /**
   * @brief Constructor with Gyper object. 
   * @details Optionally call options can be specified and a reference FASTA file.
   * 
   * @param CO Call options.
   * @param reference_fasta A reference FASTA file.
   * @return A new Gyper object.
   */
  Gyper ();
  Gyper (Options & CO);

  void add_initial_vertex();

  void add_reference_sequence_to_graph(seqan::String<seqan::Dna5> & sequence);

  // void create_reference_graph(seqan::String<char> region);

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

  int open_fasta(const char * fasta_filename);

  unsigned get_fasta_index_id(const char * id);

  void open_vcf(const char * fasta_filename);

  int read_vcf_record();

  void open_tabix(const char * fasta_filename);

  bool read_tabix_record();

  /** @brief Graph is a SeqAn partial order graph */
  TGraph graph;

  /** @brief The argument options */
  Options CO;

  /** @brief The fasta index of the reference genome */
  seqan::FaiIndex fasta_index;

  /** @brief A VCF file with all genetic variants */
  seqan::VcfFileIn vcf_file;

  /** @brief An index for BCF files, in case variants are stored in BCF instead of VCF */
  seqan::Tabix tabix_file;

  /** @brief The current VCF record. This variable is used with either VCF or BCF files. */
  seqan::VcfRecord vcf_record;


 private:
  /**
   * @brief Get the number of exons for the selected gene.
   * @details [long description]
   * @return The number of exons this HLA gene has.
   */
  unsigned get_number_of_exons(void);

  std::string get_HLA_base_path(void);

  void add_HLA_intron(void);

  void add_FASTA_region(bool add_bitstrings, int feature_number = 0, bool intron_region = false, bool p3_region = false, bool p5_region = false);
};

#endif

#ifndef __GYPER_HPP__
#define __GYPER_HPP__
// #include "graph.hpp"
#include <iostream>
#include <utility>
// #include <unordered_map>
#include <unordered_set>

#include <boost/unordered/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>

#include <gyper/options.hpp>
#include <gyper/constants.hpp>
#include <gyper/vertex_label.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>
#include <seqan/bam_io.h>
#include <seqan/stream.h>

namespace gyper
{

// Let's typesteal some things
typedef seqan::Graph<seqan::Directed<void, seqan::WithoutEdgeId> > TGraph;
typedef seqan::VertexDescriptor<TGraph>::Type TVertex;
typedef seqan::EdgeDescriptor<TGraph>::Type TEdge;

// struct KmerLabels
// {
//   TVertex start_vertex;
//   TVertex end_vertex;
//   boost::dynamic_bitset<> id_bits;
// };

/**
 * \brief Gyper's constructor class
 * \details The one class to rule them all. Each process should only have one Gyper instance.
 */

class Constructor
{
 public:
  /*********
   * BASIC *
   *********/
  TGraph graph;                        /** \brief Graph is a SeqAn partial order graph. */
  Options CO;                          /** \brief The Gyper call options. */
  seqan::GenomicRegion genomic_region; /** \brief The graph's genomic region. */
  TVertex head;                        /** \brief The newest reference vertex, i.e. the one with the highest order. */

  /**
   * \brief   Constructor for the constructor (heh). 
   * \details Optionally call options can be specified.
   * \param   CO Call options.
   * \return  A new Constructor object with default options.
   * 
   * \pre     Valid options object, if specified.
   * \post    Non-null constructor object.
   */
  Constructor ();
  Constructor (Options & CO);
  void set_genomic_region(const char * region);

  /*************
   * FASTA I/O *
   *************/
  bool fasta_set = false; 
  seqan::FaiIndex fasta_index;  /** \brief The fasta index of the reference genome */
  bool read_reference_genome(const char * fasta_filename);
  void extract_reference_sequence(void);
  void extract_reference_sequence(const char * region);


  /***********
   * VCF I/O *
   ***********/
  bool vcf_set = false; 
  seqan::VcfRecord vcf_record = seqan::VcfRecord();  /** \brief The current VCF record. This variable is used with BCF files. */
  seqan::Tabix tabix_file;                           /** \brief An index for BCF files, in case variants are stored in BCF format. */
  void open_tabix(const char * fasta_filename);
  void open_tabix(const char * tabix_filename, const char * region);
  bool read_tabix_record();
  bool read_tabix_region();
  bool read_first_tabix_region();


  /**********************
   * GRAPH CONSTRUCTION *
   **********************/
  std::vector<VertexLabel> vertex_labels;         /** \brief Vector of all nodes (vertex) labels. */
  seqan::String<seqan::Dna5> reference_sequence;  /** \brief The reference genome region we extracted from the entire thing. */
  TVertex insert_reference_vertex(unsigned const & order, seqan::String<seqan::Dna> const & value, TVertex const & previous_vertex);
  TVertex insert_reference_vertex(unsigned const & order, seqan::String<seqan::Dna> const & value);
  
  bool add_first_reference_sequence();
  void add_last_reference_sequence();
  void insert_multiple_vertexes(unsigned const & starting_pos, seqan::String<seqan::Dna> const & current_dna);
  void add_reference_sequence_preceding_a_point(unsigned const & point);
  void add_sequence_preceding_a_vcf_record(void);
  void add_vcf_record_to_graph(void);
  // TVertex add_sequence_preceding_a_vcf_record(TVertex const & current_vertex);


  /***************************
   * THE ONLY THING YOU NEED *
   ***************************/
  /**
   * \brief Given a reference, variants, and a genomic region, construct a partial order graph for that region that represents all possible haplotypes.
   * 
   * \pre   Valid location of a tabix-indexed BCF file and a fai-indexed FASTA file.
   * \post  A partial order graph which represent all possible haplotypes.
   */
  void construct_graph(const char * reference_filename, const char * vcf_filename, const char * region);
};

// Allow writing some information about the class to a stream
std::ostream &operator<<(std::ostream &os, Constructor const &c);

} // namespace gyper

#endif // __GYPER_HPP__

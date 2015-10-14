#ifndef __GRAPH_HPP_INCLUDED__
#define __GRAPH_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include <stdio.h>
#include <cstddef>
#include <string>
#include <ctime>

#include "constants.hpp"

#include <boost/unordered/unordered_set.hpp>

#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/seq_io.h>

#include <boost/unordered/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/dynamic_bitset.hpp>


using namespace seqan;


namespace seqan
{

inline
std::size_t hash_value(String<Dna> const& s)
{
  std::size_t hash_val = 0;

  for (Iterator<String<Dna> const>::Type it = begin(s) ; it != end(s) ; ++it)
  {
    hash_val = hash_val * 4 + ordValue(*it);

  }

  return hash_val;
}

inline
std::size_t hash_value(String<char> const & s)
{
  return boost::hash_range(begin(s), end(s));
}

}


struct VertexLabels {
  int level;
  Dna dna;
};

// For graph
typedef Graph<Directed<void, WithSourceId> > TGraph;
typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
typedef Size<TGraph>::Type TSize;


struct KmerLabels
{
  TVertexDescriptor start_vertex;
  TVertexDescriptor end_vertex;
  boost::dynamic_bitset<> id_bits;
};


struct GraphWithLabels {
  TGraph graph;
  std::vector<VertexLabels> vertex_labels;
};

struct ExactBacktracker {
  std::vector<TVertexDescriptor> nodes;
  std::vector<bool> match; // True if diagonal match, false otherwise
};


inline bool
operator==(const VertexLabels &lhs, const VertexLabels &rhs)
{
  return lhs.level == rhs.level && lhs.dna == rhs.dna;
}


inline bool
operator<(const VertexLabels &lhs, const VertexLabels &rhs)
{
  return lhs.level < rhs.level || (lhs.level == rhs.level && lhs.dna < rhs.dna);
}


inline std::size_t
hash_value(const VertexLabels &a)
{
  std::size_t seed = 0;
  boost::hash_combine(seed, a.level);
  boost::hash_combine(seed, ordValue(a.dna));
  return seed;
}


typedef boost::unordered_map<DnaString, std::vector<KmerLabels> > TKmerMap;
typedef boost::unordered_map<DnaString, boost::dynamic_bitset<> > TKmerMapSimple;

struct callOptions
{
 public:
  CharString beta_list;
  std::vector<double> beta;
  int bpQclip;
  int bpQskip;
  int number_of_exons;
  int kmer;
  CharString gene;
  CharString minSeqLen_list;
  CharString outputFolder;
  CharString vcfOutputFolder;
  CharString bamFile;
  CharString bam2;
  CharString bam3;
  CharString bam4;
  std::vector<unsigned> minSeqLen;
  CharString vcfFile;
  
  bool verbose;
  bool align_all_reads;
  bool thousand_genomes;
  int read_gap;

  void output(){
    std::cout << "CO.bamFile " << bamFile << std::endl;
    std::cout << "CO.bam2 " << bam2 << std::endl;
    std::cout << "CO.bam3 " << bam3 << std::endl;
    std::cout << "CO.bam4 " << bam4 << std::endl;
    std::cout << "CO.beta_list " << beta_list << std::endl;
    std::cout << "CO.beta.size() " << beta.size() << std::endl;
    std::cout << "CO.bpQclip " << bpQclip << std::endl;
    std::cout << "CO.bpQskip " << bpQskip << std::endl;
    std::cout << "CO.gene " << gene << std::endl;
    std::cout << "CO.number_of_exons " << number_of_exons << std::endl;
    std::cout << "CO.kmer " << kmer << std::endl;
    std::cout << "CO.minSeqLen_list " << minSeqLen_list << std::endl;
    std::cout << "CO.minSeqLen.size() " << minSeqLen.size() << std::endl;
    std::cout << "CO.outputFolder " << outputFolder << std::endl;
    std::cout << "CO.vcfFile " << vcfFile << std::endl;
    std::cout << "CO.vcfOutputFolder" << vcfOutputFolder << std::endl;
    std::cout << "CO.verbose " << verbose << std::endl;
    std::cout << "CO.align_all_reads " << align_all_reads << std::endl;
    std::cout << "CO.thousand_genomes " << thousand_genomes << std::endl;
    std::cout << "CO.read_gap " << read_gap << std::endl;
  };

callOptions():
  beta_list("0.6"), bpQclip(30), bpQskip(25), number_of_exons(4), kmer(15), gene("DQA1"), minSeqLen_list("60"), outputFolder(),
  vcfOutputFolder(), verbose(false), align_all_reads(false), thousand_genomes(false), read_gap(1000) {}
};


void
align_sequence (DnaString & my_sequence,
                boost::dynamic_bitset<> & qual,
                TGraph const & graph,
                std::vector<VertexLabels> & vertex_vector,
                String<TVertexDescriptor> & order,
                std::vector<ExactBacktracker> & backtracker,
                std::vector<ExactBacktracker> & reverse_backtracker,
                boost::unordered_set<TVertexDescriptor> const & free_nodes,
                std::vector<TVertexDescriptor> & matching_vertices,
                std::vector<TVertexDescriptor> & reverse_matching_vertices
               );

boost::dynamic_bitset<>
align_sequence_kmer (String<Dna> & my_sequence,
                     String<char> & qual,
                     unsigned const & id_numbers,
                     TKmerMap & kmer_map,
                     std::vector<VertexLabels> & vertex_vector,
                     int const & kmer_size
                    );

CharString myExtractTagValue(String<char> &tags);

void
createGenericGraph(callOptions & CO,
                   TGraph & graph,
                   std::vector<VertexLabels> & vertex_vector,
                   std::vector<std::string> & ids,
                   boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                   boost::unordered_set<TVertexDescriptor> & free_nodes,
                   String<TVertexDescriptor> & order
                  );

void
createDqa1Graph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createDqb1Graph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createDrb1Graph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createHlaaGraph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createHlabGraph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);

void
createHlacGraph(TGraph & graph,
                std::vector<VertexLabels> & vertex_vector,
                std::vector<std::string> & ids,
                boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > & edge_ids,
                boost::unordered_set<TVertexDescriptor> & free_nodes,
                String<TVertexDescriptor> & order);


#include "graph_builder.hpp"
#include "graph_align.hpp"
#include "graph_io.hpp"
#include "graph_kmerify.hpp"

#endif

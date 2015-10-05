#ifndef __GRAPH_HPP_INCLUDED__
#define __GRAPH_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include <stdio.h>
#include <ctime>

#include "constants.hpp"
#include "graph_align.hpp"
#include "graph_builder.hpp"
#include "graph_io.hpp"

#include <boost/unordered/unordered_set.hpp>

#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include <seqan/seq_io.h>

using namespace seqan;


struct callOptions
{
 public:
  CharString beta_list;
  std::vector<double> beta;
  int bpQclip;
  int bpQskip;
  CharString bamFile;
  CharString bam2;
  CharString bam3;
  CharString bam4;
  CharString gene;
  int number_of_exons;
  CharString minSeqLen_list;
  std::vector<unsigned> minSeqLen;
  // std::vector<int> minSeqs;
  CharString outputFolder;
  CharString vcfFile;
  CharString vcfOutputFolder;
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
  beta_list("0.6"), bpQclip(30), bpQskip(25), number_of_exons(4), gene("DQA1"), minSeqLen_list("60"), outputFolder(), vcfOutputFolder(), verbose(false), align_all_reads(false), thousand_genomes(false), read_gap(1000) {}
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

#endif

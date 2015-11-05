#ifndef __GRAPH_IO_HPP_INCLUDED__
#define __GRAPH_IO_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <time.h>

#include "graph.hpp"

#include <seqan/vcf_io.h>
#include <seqan/sequence.h> 

void printIdMap(boost::unordered_map<std::string, long> &my_id_map);
void writeDotFile(TGraph graph);
void printGraph(TGraph graph);

void
writeVcfFile(seqan::CharString & file_name,
             std::string const & pn,
             std::vector<std::string> ids,
             std::vector<std::vector<double> > & seq_scores,
             double const & max_score
            );

void
writeVcfFile(seqan::CharString & file_name,
             std::string const & pn,
             std::vector<std::string> ids,
             std::vector<std::vector<double> > & seq_scores,
             double const & max_score,
             unsigned const & minSeqLen,
             double const & beta,
             int const & bpQclip,
             int const & bpQskip
            );

void
testVcfFile(seqan::CharString file_name,
            std::string pn,
            std::vector<std::string> ids,
            std::vector<std::vector<double> > seq_scores,
            double max_score
            );

#endif

#ifndef __GRAPH_IO_HPP_INCLUDED__
#define __GRAPH_IO_HPP_INCLUDED__

#define SEQAN_NO_GLOBAL_EXCEPTION_HANDLER

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <time.h>

#include "graph_align.hpp"
#include "graph_builder.hpp" // Needs the templates defined there

#include <seqan/vcf_io.h>
#include <seqan/sequence.h> 

void printIdMap(boost::unordered_map<std::string, long> &my_id_map);
void writeDotFile(TGraph graph);
void printGraph(TGraph graph);

void
writeVcfFile(CharString & file_name,
             std::string const & pn,
             std::vector<std::string> ids,
             std::vector<std::vector<double> > & seq_scores,
             double const & max_score
            );

void
writeVcfFile(CharString & file_name,
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
testVcfFile(CharString file_name,
            std::string pn,
            std::vector<std::string> ids,
            std::vector<std::vector<double> > seq_scores,
            double max_score
            );

#endif

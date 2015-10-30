#ifndef GYPER_OPTIONS_HPP
#define GYPER_OPTIONS_HPP

#include <seqan/basic.h>
#include <seqan/arg_parse.h>

class Options
{
 public:
  double beta;                         /** @brief A parameter which scales heterozygous scores down to better adjust for zygosity. */
  int k;                               /** @brief Length of a k-mer. */
  int min_kmers;
  seqan::String<char> alignment_file;
  seqan::String<char> gene;
  seqan::String<char> output_folder;
  bool verbose;
  bool exon_2_and_3;
  bool vcf;
  bool bias_check;
  bool thousand_genomes;
  int read_gap;

  void print_options()
  {
    std::cout << "Option beta = " << beta << std::endl;
    std::cout << "Option alignment_file = " << alignment_file << std::endl;
    std::cout << "Option gene = " << gene << std::endl;
    std::cout << "Option k = " << k << std::endl;
    std::cout << "Option min_kmers = " << min_kmers << std::endl;
    std::cout << "Option verbose = " << verbose << std::endl;
    std::cout << "Option exon_2_and_3 = " << exon_2_and_3 << std::endl;
    std::cout << "Option vcf = " << vcf << std::endl;
    std::cout << "Option bias_check = " << bias_check << std::endl;
    std::cout << "Option thousand_genomes = " << thousand_genomes << std::endl;
    std::cout << "Option read_gap = " << read_gap << std::endl;
  }

  Options()
  {
    beta = 0.6;
    k = 15;
    min_kmers = 1;
    gene = "DQA1";
    verbose = false;
    exon_2_and_3 = false;
    vcf = false;
    bias_check = false;
    thousand_genomes = false;
    read_gap = 1000;
  }
};

/**
 * @brief Parses command line options
 * @details [long description]
 * 
 * @param CO A option object which will be inserted data into
 * @param parser [description]
 * @param argc [description]
 * @param argv [description]
 * @return A value indicating if parsing the command line was successful or not.
 */
seqan::ArgumentParser::ParseResult parse_command_line_options(Options & CO, seqan::ArgumentParser & parser, int argc, char const ** argv);
#endif

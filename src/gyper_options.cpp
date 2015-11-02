#include "gyper_options.hpp"

seqan::ArgumentParser::ParseResult parse_command_line_options(Options & CO, seqan::ArgumentParser & parser, int argc, char const ** argv)
{
  using namespace seqan;

  setShortDescription(parser, "Genotypes a given BAM file.");
  setVersion(parser, "0.1");
  setDate(parser, "November 2015");

  // Main options
  addSection(parser, "Main options");
  addOption(parser, ArgParseOption("k", "k", "The length of the k-mers.", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("mk", "min_kmers", "The minimum amount of k-mers.", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("v", "verbose", "Verbosity flag"));
  addOption(parser, ArgParseOption("23", "exon_2_and_3", "Only use graph with exons 2 and 3."));
  addOption(parser, ArgParseOption("c", "vcf", "Print VCF outut."));
  addOption(parser, ArgParseOption("bc", "bias_check", "If used Gyper will use every read pair in the BAM file provided. Should only be used for very small BAMs, unless you have huge amount of RAM and patience."));
  addOption(parser, ArgParseOption("1k", "thousand_genomes", "Use reference from the 1000 Genomes project."));
  addOption(parser, ArgParseOption("b", "beta", "A scoring parameter for heterozygous scores. Higher beta means more likely to pick heterozugous solutions.", ArgParseArgument::INTEGER, "INTEGER"));
  setDefaultValue(parser, "beta", 0.6);
  addOption(parser, ArgParseOption("rg", "read_gap", "The amount af bp to extend the read gap (i.e. the gap where reads are considered.) ", ArgParseArgument::INTEGER, "INT"));

  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // Only extract options if the program will continue after parseCommandLine()
  if (res != ArgumentParser::PARSE_OK)
  {
    return res;
  }

  getOptionValue(CO.k, parser, "k");
  getOptionValue(CO.min_kmers, parser, "min_kmers");
  getOptionValue(CO.beta, parser, "beta");
  getOptionValue(CO.bias_check, parser, "bias_check");
  getOptionValue(CO.thousand_genomes, parser, "thousand_genomes");
  getOptionValue(CO.read_gap, parser, "read_gap");

  // First parameter is the gene to genotype, second is a bam file (or a file containing a list of bam files)
  getArgumentValue(CO.gene, parser, 0);
  getArgumentValue(CO.alignment_file, parser, 1);
  CO.verbose = isSet(parser, "verbose");
  CO.exon_2_and_3 = isSet(parser, "exon_2_and_3");
  CO.vcf = isSet(parser, "vcf");

  // Remove the HLA- part from genes if it is there
  if (CO.gene == "HLA-A")
  {
    CO.gene = "HLAA";
  }
  else if (CO.gene == "HLA-B")
  {
    CO.gene = "HLAB";
  }
  else if (CO.gene == "HLA-C")
  {
    CO.gene = "HLAC";
  }
  else if (CO.gene == "HLA-DQA1")
  {
    CO.gene = "DQA1";
  }
  else if (CO.gene == "HLA-DQB1")
  {
    CO.gene = "DQB1";
  }
  else if (CO.gene == "HLA-DRB1")
  {
    CO.gene = "DRB1";
  }

  if (CO.verbose)
  {
    CO.print_options();
  }

  return ArgumentParser::PARSE_OK;
}

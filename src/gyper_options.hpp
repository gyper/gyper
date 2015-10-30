#include <seqan/basic.h>

class Options
{
 private:
  /*
   * Instance variables
   **/

  unsigned number_of_exons;
  seqan::String<char> short_gene_name;

 public:
  double beta;
  int k;
  int min_kmers;
  seqan::String<char> gene;
  seqan::String<char> output_folder;
  seqan::String<char> alignment_file;
  unsigned minSeqLen;
  
  bool verbose;
  bool exon_2_and_3;
  bool vcf;
  bool bias_check;
  bool thousand_genomes;
  int read_gap;

  void print_options()
  {
    std::cout << "Option bam_file = " << bamFile << std::endl;
    std::cout << "Option gene = " << gene << std::endl;
    std::cout << "Option k = " << k << std::endl;
    std::cout << "Option min_kmers = " << min_kmers << std::endl;
    std::cout << "Option minSeqLen_list = " << minSeqLen_list << std::endl;
    std::cout << "Option minSeqLen.size() = " << minSeqLen.size() << std::endl;
    std::cout << "Option outputFolder = " << outputFolder << std::endl;
    std::cout << "Option vcfFile = " << vcfFile << std::endl;
    std::cout << "Option vcfOutputFolder = " << vcfOutputFolder << std::endl;
    std::cout << "Option verbose = " << verbose << std::endl;
    std::cout << "Option exon_2_and_3 = " << exon_2_and_3 << std::endl;
    std::cout << "Option vcf = " << vcf << std::endl;
    std::cout << "Option bias_check = " << bias_check << std::endl;
    std::cout << "Option thousand_genomes = " << thousand_genomes << std::endl;
    std::cout << "Option read_gap = " << read_gap << std::endl;
  }

  void print_all()
  {
    std::cout << "number_of_exons = " << number_of_exons << std::endl;
  }

  callOptions()
  {
    beta_list = "0.6";
    number_of_exons(4), kmer(15), min_kmers(1), gene("DQA1"), minSeqLen_list("60"), outputFolder(),
  vcfOutputFolder(), verbose(false), exon_2_and_3(false), vcf(false), bias_check(false), thousand_genomes(false), read_gap(1000) {}
  }

  void set_gene(seqan::String<char> new_gene)
  {
    switch (new_gene)
    {
      case "HLAA":
      case "HLA-A":
        number_of_exons = 8;
      case "HLAB":
      case "HLA-B":

      case "HLAC":
      case "HLA-C":


      case "DQA1":
      case "HLA-DQA1":
        number_of_exons = 4;
        break;
      case "DQB1":
      case "HLA-DQB1":

      case "DRB1":
      case "HLA-DRB1":

    }
  }
};


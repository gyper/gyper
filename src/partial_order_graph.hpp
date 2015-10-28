#include "graph.hpp"

class POGraph
{
  private:
    /*
      Private instance variables
    */
    callOptions CO;
    std::vector<VertexLabels> vertex_vector;
    std::vector<std::string> ids;
    boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;
    boost::unordered_set<TVertexDescriptor> free_nodes;
    String<TVertexDescriptor> order;

    // Testing variables
    int width;
    int height;

    /*
     Private functions
    */

    // Pre conditions : CO.gene is some HLA gene
    // Post conditions: A reference graph with all known HLA alles for gene in CO.gene 
    void create_HLA_exon_2_and_3_graph(callOptions & CO);


  public:
    POGraph ();
    POGraph (callOptions & CO);

    /*
     Public variables
    */
    TGraph graph; // graph is a SeqAn graph

    void set_values (int, int);
    int area ();
};

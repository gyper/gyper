#include "graph.hpp"

/**
 * @brief Gyper's main class
 * @details [long description]
 */

class Gyper
{
 private:
  std::vector<VertexLabels> vertex_vector;
  std::vector<std::string> ids;
  boost::unordered_map< std::pair<TVertexDescriptor, TVertexDescriptor>, boost::dynamic_bitset<> > edge_ids;
  boost::unordered_set<TVertexDescriptor> free_nodes;
  String<TVertexDescriptor> order;

  // Testing variables
  int width;
  int height;

  /**
   * @brief Creates a single HLA graph
   * @details Creates an partial order graph which represents a HLA gene specified by the call option
   */
  void create_HLA_graph();

  public:
    Gyper ();
    Gyper (callOptions & CO);

    /*
     Public variables
    */
    TGraph graph; // graph is a SeqAn graph
    callOptions CO; // The argument options

    void set_values (int, int);
    int area ();
};

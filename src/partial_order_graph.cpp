#include "partial_order_graph.hpp"

/*
 Constructors
*/

// Post conditions: Empty partial order graph
POGraph::POGraph ()
{
	TGraph graph;
  POGraph::create_HLA_exon_2_and_3_graph(CO);
}

POGraph::POGraph (callOptions & CO)
{
  // Constructor with specific callOptions
  POGraph::POGraph::CO = CO;
  POGraph();
}

void
POGraph::create_HLA_exon_2_and_3_graph(callOptions & CO)
{
  int number_of_exons = CO.number_of_exons;

  std::stringstream base_path;
  base_path << gyper_SOURCE_DIRECTORY << "/data/haplotypes/hla/references/" << CO.gene << "/";
  TVertexDescriptor begin_vertex;

  {
    std::string p3_string = base_path.str();
    p3_string.append("p3.fa");

    if (CO.verbose)
    {
      std::cout << "Adding utr    " << p3_string << std::endl;
    }

    const char* alignment_file_p3 = p3_string.c_str();
    graph = createGraph(alignment_file_p3, vertex_vector, ids, begin_vertex);
  }
  
  TVertexDescriptor new_begin_vertex;

  std::string tmp_string;
  std::string extension = ".fa";

  // First exon
  {
    tmp_string = base_path.str();
    std::string exon = "e";
    std::string feature_number = std::to_string(number_of_exons);
    exon.append(feature_number);
    exon.append(extension);
    tmp_string.append(exon);

    // if (CO.verbose)
    // {
    //   std::cout << "Adding exon   " << tmp_string << std::endl;
    // }

    const char* alignment_file = tmp_string.c_str();
    free_nodes.insert(begin_vertex);

    if (CO.verbose)
    {
      std::cout << "Adding exon   " << tmp_string << " as intron" << std::endl;
    }

    extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);

    --number_of_exons;
  }

  while (number_of_exons >= 1)
  {
    // Intron
    {
      tmp_string = base_path.str();
      std::string intron = "i";
      std::string feature_number = std::to_string(number_of_exons);
      intron.append(feature_number);
      intron.append(extension);
      tmp_string.append(intron);

      if (CO.verbose)
      {
        std::cout << "Adding intron " << tmp_string << std::endl;
      }

      const char* alignment_file = tmp_string.c_str();
      free_nodes.insert(begin_vertex);
      extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);
    }

    // Exon
    {
      tmp_string = base_path.str();
      std::string exon = "e";
      std::string feature_number = std::to_string(number_of_exons);
      exon.append(feature_number);
      exon.append(extension);
      tmp_string.append(exon);

      // if (CO.verbose)
      // {
      //   std::cout << "Adding exon   " << tmp_string << std::endl;
      // }

      const char* alignment_file = tmp_string.c_str();
      free_nodes.insert(begin_vertex);
      
      if (number_of_exons == 2 || number_of_exons == 3 || number_of_exons == 4)
      {
        if (CO.verbose)
        {
          std::cout << "Adding exon   " << tmp_string << std::endl;
        }

        extendGraph(graph, alignment_file, vertex_vector, edge_ids, new_begin_vertex, begin_vertex);
      }
      else
      {
        if (CO.verbose)
        {
          std::cout << "Adding exon   " << tmp_string << " as intron" << std::endl;
        }

        extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);
      }
    }

    --number_of_exons;
  }

  {
    // Final UTR
    tmp_string = base_path.str();
    std::string utr = "5p";
    utr.append(extension);
    tmp_string.append(utr);

    if (CO.verbose)
    {
      std::cout << "Adding utr    " << tmp_string << std::endl;
    }

    const char* alignment_file = tmp_string.c_str();
    free_nodes.insert(begin_vertex);
    extendGraph(graph, alignment_file, vertex_vector, new_begin_vertex, begin_vertex);
  }

  topologicalSort(order, graph);
}

void POGraph::set_values (int x, int y)
{
  width = x;
  height = y;
}

int POGraph::area ()
{
  return POGraph::width*POGraph::height;
}


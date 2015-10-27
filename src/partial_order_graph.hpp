#include "graph.hpp"

class POGraph
{
	private:
		TGraph graph;
    int width;
    int height;
  public:
    void set_values (int,int);
    int area(); // {return width*height;}
};

#include "fmmtl/tree/KDBTree.hpp"
#include "fmmtl/numeric/Vec.hpp"

#include <iostream>
#include <string>

#define DIM 3

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NUM_POINTS\n";
    exit(1);
  }

  int N = atoi(argv[1]);

  typedef Vec<DIM,double> point_type;

  std::vector<point_type> points(N);
  for (int k = 0; k < N; ++k)
    points[k] = fmmtl::random<point_type>::get();

  // if (argc == 4) {
  //   int nr = atoi(argv[2]);
  //   int np = atoi(argv[3]);
  //   NDBTree<DIM> tree(points.begin(), points.end(), nr, np);
  // }
  // else {
    NDBTree<DIM> tree(points.begin(), points.end());
    auto i = 0;
    for (auto &&p : points) {
      assert(tree.query(p));
      ++i;
    }
    if (argc == 3) {
      if (atoi(argv[2]))
        tree.print_graph(1); // print nodes
      else
        tree.print_graph(0); // print edges
    }
  
}

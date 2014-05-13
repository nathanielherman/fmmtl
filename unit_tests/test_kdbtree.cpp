#include "fmmtl/tree/KDBTree.hpp"
#include "fmmtl/numeric/Vec.hpp"

#include <iostream>
#include <string>


#define DIM 2

typedef Vec<DIM,double> point_type;

void test_times(int N) {
  for (int i = 0; i < N; i += 100000) {
    std::vector<point_type> points(i);
    for (int k = 0; k < i; ++k)
      points[k] = fmmtl::random<point_type>::get();

    clock_t start = clock();
    NDBTree<DIM> tree(points.begin(), points.end());
    clock_t end = clock();
    std::cout << i << ", " << (end - start) / (CLOCKS_PER_SEC/1000) << std::endl;
  }
}

void test_interesting(int N, int argc, char **argv) {

  // std::vector<point_type> points(N);

  // srand(1);
  // for (int k = 0; k < N-100; ++k) {
  //   point_type p;
  //   for (int i = 0; i < DIM; ++i) {
  //     p[i] = (double) rand() / RAND_MAX;
  //   }
  //   points[k] = p;
  // }
  // for (int k = 0; k < 100; ++k) {
  //   points[N-1-k] = point_type(1.5+k*.01,1.5,1.5);  
  // }

  //   NDBTree<DIM> tree(points.begin(), points.end());

  //   if (argc == 3) {
  //     if (atoi(argv[2]))
  //       tree.print_graph(1); // print nodes
  //     else
  //       tree.print_graph(0); // print edges
  //   }
}

void test_normal(int N, int argc, char **argv) {

  std::vector<point_type> points(N);
  for (int k = 0; k < N; ++k)
    points[k] = fmmtl::random<point_type>::get();

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

void test_query(int N) {
  for (int i = 4000000; i < 5100000; i += 100000) {
    std::vector<point_type> points(i);
    for (int k = 0; k < i; ++k)
      points[k] = fmmtl::random<point_type>::get();
      
    NDBTree<DIM> tree(points.begin(), points.end());
      
    clock_t start = clock();
    for (auto &&p : points) {
      assert(tree.query(p));
    }
    clock_t end = clock();
    std::cout << i << ", " << (double) ((end - start) / (CLOCKS_PER_SEC/1000)) / i << std::endl;
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NUM_POINTS\n";
    exit(1);
  }

  int N = atoi(argv[1]);

  // test_normal(N, argc, argv);
  test_query(N);
  return 0;  
}

#include "fmmtl/tree/KDBTree.hpp"
#include "fmmtl/numeric/Vec.hpp"

#include <iostream>
#include <string>


#define DIM 2

enum {
  Times,
  Interesting,
  Normal,
  Query,
  QueryRange
};

#define MODE Normal

typedef Vec<DIM,double> point_type;

void test_times(int N) {
  for (int i = 100000; i < N; i += 100000) {
    std::vector<point_type> points(i);
    for (int k = 0; k < i; ++k)
      points[k] = fmmtl::random<point_type>::get();

    clock_t start = clock();
    KDBTree<DIM> tree(points.begin(), points.end());
    clock_t end = clock();
    std::cout << i << ", " << (end - start) / (CLOCKS_PER_SEC/1000) << std::endl;
  }
}

// Implements a square pattern and some outliers
void test_interesting(int N, int argc, char **argv) {
  std::vector<point_type> points(N);

  srand(1);
  for (int k = 0; k < N-100; ++k) {
    point_type p;
    for (int i = 0; i < DIM; ++i) {
      p[i] = (double) rand() / RAND_MAX;
      p[i]*=p[i];
    }
    points[k] = p;
  }
  for (int k = 0; k < 100; ++k) {
    point_type p;
    for (int i = 0; i < DIM; ++i) {
      p[i] = 1.5;
    }
    p[0] += k*.01;
    points[N-1-k] = p;
  }
  KDBTree<DIM> tree(points.begin(), points.end());

  if (argc == 3) {
    if (atoi(argv[2]))
      tree.print_graph(1); // print nodes
    else
      tree.print_graph(0); // print edges
  }
}

void test_normal(int N, int argc, char **argv) {

  std::vector<point_type> points(N);
  for (int k = 0; k < N; ++k)
    points[k] = fmmtl::random<point_type>::get();

  KDBTree<DIM> tree(points.begin(), points.end());

  // print nodes and edges
  if (argc == 3) {
    if (atoi(argv[2]))
      tree.print_graph(1); // print nodes
    else
      tree.print_graph(0); // print edges
  }
}

void test_query(int N) {
  for (int i = 100000; i < N; i += 100000) {
    std::vector<point_type> points(i);
    for (int k = 0; k < i; ++k)
      points[k] = fmmtl::random<point_type>::get();
      
    KDBTree<DIM> tree(points.begin(), points.end());
      
    clock_t start = clock();
    for (auto &&p : points) {
      assert(tree.query(p));
    }
    clock_t end = clock();
    std::cout << i << ", " << (double) ((end - start) / (CLOCKS_PER_SEC/1000)) / i << std::endl;
  }
}

void test_query_range(int N) {
  for (int i = 100000; i < N; i += 100000) {

    std::vector<point_type> points(i);
    for (int k = 0; k < i; ++k)
      points[k] = fmmtl::random<point_type>::get();
      
    KDBTree<DIM> tree(points.begin(), points.end());
    
    clock_t start = clock();
    for (int i = 0; i < 10000; ++i) {
      point_type min(.1*(i%10),.1*(i%10));
      point_type max = min;
      max[0] += 100.0/i;
      auto vec = tree.query_range(BoundingBox<DIM>(min,max));
    }
    clock_t end = clock();
    std::cout << i << ", " << (double) ((end - start) / (CLOCKS_PER_SEC/1000)) << std::endl;    
  }
}


int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NUM_POINTS\n";
    exit(1);
  }

  int N = atoi(argv[1]);

  switch(MODE) {
    case Normal:
      test_normal(N, argc, argv);
      break;
    case Times:
      test_times(N);
      break;
    case Interesting:
      test_interesting(N, argc, argv);
      break;
    case Query:
      test_query(N);
      break;
    case QueryRange:
      test_query_range(N);
      break;
  }

  return 0;  
}

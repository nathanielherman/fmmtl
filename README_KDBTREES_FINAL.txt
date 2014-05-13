KDBTree Implemenation
=============================

The implementation of the KDBTree is based on this paper http://dl.acm.org/citation.cfm?id=582321 and can be found under ffmtl/tree/KDBTree.hpp

Our implementation has been tested and works for 2D and 3D points.     

The KDBTree can be instantiated by:

KDBTree<DIM> tree(POINT_RANGE_BEGIN, POINT_RANGE_END);
KDBTree<DIM> tree(POINT_RANGE_BEGIN, POINT_RANGE_END, n_crit_region, n_crit_point);

In the first case we pass in a range of points of DIM dimension. 

In the second case we also specify the max number of regions we can have in a region page before splitting and the max number of points we can have in a point page.

If those are not specified we choose them so that a region page and a point page can fit in cache.

### Run - Works for all modes
./test_kdbtree NUM_POINTS

### Visualization - Works for 2D and 3D in Normal Mode, Interesting Mode

./test_kdbtree NUM_POINTS 1
prints a node file of the tree
./test_kdbtree NUM_POINTS 0
prints the edge file of the tree

These files can be fed into the CS207 visualizer to produce a cool visualization of the results as also seen in our report.

### Dimension

Inside the test_kdbtree.hpp we specify DIM at the top. 
If we want to print the nodes and the edges of the bounding boxes DIM has to be 2 or 3.
Other than the visualization the KDBTree should work for an arbitrary number of dimensions.

### MODES

Times: tests loading times for different tree sizes
Interesting: generates points from a uniform squared distribution and adds some outliers
Normal: creates uniformly random points between (0) and (1)
Query: tests the speed of single point queries
QueryRange: tests the speed of 100 point range queries in different tree sizes

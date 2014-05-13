#pragma once
/** @file NDTree
 * @brief General class representing a {1D,2D,3D,4D}-BTree.
 */

#include <vector>
#include <bitset>
#include <algorithm>

#include <iostream>
#include <iomanip>

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/permutation_iterator.hpp>
using boost::iterator_adaptor;

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/numeric/Vec.hpp"
#include "BoundingBox.hpp"

#define CACHE_SZ 4096

//! Class for tree structure
template <unsigned DIM>
struct NDBTree {


  // type declarations
  typedef Vec<DIM,double> point_type;
  typedef BoundingBox<DIM> bounding_box_type;

  struct RegionPage;
  struct Page {
    bool isRegionPage;
    uint8_t splittingDomain;
    uint16_t pidx;
    RegionPage *parent;
    
    Page() : isRegionPage(false), splittingDomain(0), pidx(0), parent(NULL) {}
    
    virtual ~Page() {}
  };
  
  struct region_type {
    bounding_box_type box;
    Page *page;
  };

  struct RegionPage : public Page {
    std::vector<region_type> children;
    RegionPage() : Page(), children() { this->isRegionPage = true; }
  };

  struct PointPage : public Page {
    std::vector<point_type> points;
    PointPage() : Page(), points() { this->isRegionPage = false; }
  };

  // Precomputed n_crit_region_ and n_crit_point_
  static constexpr unsigned n_crit_region = CACHE_SZ / (sizeof(region_type)) - CACHE_SZ/100;
  static constexpr unsigned n_crit_point = CACHE_SZ / (sizeof(point_type)) - CACHE_SZ/100;
  
  unsigned n_crit_region_;
  unsigned n_crit_point_;
  
  template <typename PointIter>
  /** NDBTree constructor
   * @param[in] first, last Insert points in range [first, last)
   * @param[in] max_reg Optional, maximum regions per region page
   * @param[in] max_pt Optional, maximum points per point page
   * @pre No duplicate points in [first, last)
   */
  NDBTree(PointIter first, PointIter last, unsigned max_reg = n_crit_region, unsigned max_pt = n_crit_point)
      : n_crit_region_(max_reg), n_crit_point_(max_pt) {
    std::cout << "n_crit_region " << n_crit_region_ << ", n_crit_point " << n_crit_point_ << std::endl;
    
    // create empty root page
    PointPage *pp = new PointPage();
    pointPages_.push_back(pp);
    root = pp;
    rootBox = get_boundingbox(first, last);
    
    insert_range(first, last);
  }

  /** NDBTree destructor */
  ~NDBTree() {
    for (PointPage *p : pointPages_) {
      delete p;
    }
    for (RegionPage *p : regionPages_) {
      delete p;
    }
  }

  /************************************
   ********** PRINTING HELPERS ********
   ************************************/
  
   void print(Page *p) {
      if (!p->isRegionPage) {
          PointPage *pp = dynamic_cast<PointPage*> (p);
          std::cout << "Point Page with bounding box" << std::endl;
          for (auto&& point: pp->points) {
   // std::cout << point << std::endl;
          }
   // std::cout << std::endl;
        return;
      }
      RegionPage *rp = dynamic_cast<RegionPage*> (p);
      for (auto&& region: rp->children) {
        std::cout << "Bounding box: " << region.box << std::endl;
        print(region.page);
      }
  }

  void print() {
    // std::cout << "\n\n\n\n Printing Tree " << std::endl;
    // std::cout << "Root bounding box" << rootBox << std::endl;
    print(root);
    // std::cout << "End of tree\n\n" << std::endl;
  }

  // @brief helper for print_graph
  void print_edge(size_t p1, size_t p2) {
    std::cout << p1 << " " << p2 << std::endl;
  }
  // @brief helper for print_graph
  void print_point_3d(point_type p) {
    std::cout << p << std::endl; 
  }

  /* prints bounding box
   * @pre DIM is 2 or 3
   * @post if DIM = 2, old_num_nodes + 4 = new_num_nodes
   * @post if DIM = 3, old_num_nodes + 8 = new_num_nodes
   * @param nodes tells us whether to print nodes or edges
   */

  size_t print_box(bounding_box_type box, size_t num_nodes, bool nodes) {
    if (DIM == 2) {
      point_type p1 = box.min();
      point_type p2 = box.min();
      p1[0] = box.max()[0];
      p2[1] = box.max()[1];

      size_t min_n = num_nodes++;
      size_t max_n = num_nodes++;
      size_t p1_n = num_nodes++;
      size_t p2_n = num_nodes++;

      if (nodes) {
        std::cout << box.min() << " 0 " << std::endl;
        std::cout << box.max() << " 0 " <<std::endl;
        std::cout << p1 << " 0 " << std::endl;
        std::cout << p2 << " 0 " << std::endl;
      }
      else {
        print_edge(min_n, p1_n); print_edge(min_n, p2_n);
        print_edge(max_n, p1_n); print_edge(max_n, p2_n);
      } 
    }
    else if (DIM == 3) {
      point_type p1 = box.min(); point_type p2 = box.min(); point_type p3 = box.min();
      point_type p4 = box.min(); point_type p5 = box.min(); point_type p6 = box.min();

      p1[0] = box.max()[0];
      p2[0] = box.max()[0]; p2[1] = box.max()[1];
      p3[0] = box.max()[0]; p3[2] = box.max()[2];
      p4[1] = box.max()[1]; p4[2] = box.max()[2];
      p5[1] = box.max()[1];
      p6[2] = box.max()[2];

      size_t min_n = num_nodes++;
      size_t max_n = num_nodes++;
      size_t p1_n = num_nodes++; size_t p2_n = num_nodes++; size_t p3_n = num_nodes++;
      size_t p4_n = num_nodes++; size_t p5_n = num_nodes++; size_t p6_n = num_nodes++;

      if (nodes) {
        print_point_3d(box.min()); print_point_3d(box.max());
        print_point_3d(p1); print_point_3d(p2); print_point_3d(p3);
        print_point_3d(p4); print_point_3d(p5); print_point_3d(p6);
      }
      else {
        print_edge(min_n, p1_n); print_edge(min_n, p5_n); print_edge(min_n, p6_n);
        print_edge(max_n, p2_n); print_edge(max_n, p3_n); print_edge(max_n, p4_n);
        print_edge(p1_n, p2_n); print_edge(p1_n, p3_n); print_edge(p2_n, p5_n);
        print_edge(p5_n, p4_n); print_edge(p3_n, p6_n); print_edge(p4_n, p6_n);
      } 
    }
    return num_nodes;
  }

  size_t print_graph_rest(Page *p, size_t num_nodes, bool nodes) {
    if (p->isRegionPage) {
      RegionPage *rp = dynamic_cast<RegionPage*> (p);
      for (auto&& region: rp->children) {
        num_nodes = print_box(region.box, num_nodes, nodes);
        num_nodes = print_graph_rest(region.page, num_nodes, nodes);
      }
    }
    return num_nodes;
  }

  /* prints all the points and adds a zero coord at the end so that it
   * works with the visualizer
   */
  void print_points() {
    for (auto &&page : pointPages_) {
      for (auto && point : page->points) {
        assert(page->points.size() <= n_crit_point_);
        std::cout << point << " 0" << std::endl;
      }
    }
  }

  // prints nodes and edges for the visualizer
  void print_graph(bool nodes) {
    if (DIM == 2 || DIM == 3) {
      size_t num_nodes = 0;
      num_nodes = print_box(rootBox, num_nodes, nodes);
      print_graph_rest(root, num_nodes, nodes);
      if (nodes) {
        print_points();
      }
    } 
  }
  
  private:

  /** Inserts the given range of points into the tree */
  template <typename PointIter>
  void insert_range(PointIter p_first, PointIter p_last) {
    for (auto it = p_first; it != p_last; ++it) {
      if (!insert(*it)) {
        // inserted duplicate point // TODO: notify user
        assert(0);
      }
    }
  }
  
  /** Insert a point into the tree
   * @param[in] p the point to add
   * @return false if p is already in the tree, true otherwise
   * @post query(p)==true
   * Complexity: O(logn)
   */
  bool insert(point_type p) {
    // we create the root in the constructor
    assert(root != NULL);
    assert(rootBox.contains(p));
    
    PointPage &page = *_query(p);
    // this point was already inserted
    if (inPointPage(p, &page)) {
      return false;
    }

    page.points.push_back(p);
    
    if (page.points.size() <= n_crit_point_) {
      return true;
    }
    
    // need to split this point page :(
    split(page);
    
    return true;
  }

  /** Returns where to split the given region page 
   * (using p.splittingDomain as the dimension to split on) 
   */
  double calcRegionSplit(RegionPage& p) {
    auto& regions = p.children;
    auto dim = p.splittingDomain;
    //std::nth_element(regions.begin(), median, regions.end(), 
    std::sort(regions.begin(), regions.end(),
	      [=] (const region_type r1, const region_type r2) { 
		return r1.box.min()[dim] < r2.box.min()[dim]; 
	      });
    auto median = regions.begin() + regions.size()/2;
    auto minPt = regions[0].box.min()[p.splittingDomain];
    while (median != regions.end() 
	   && (*median).box.min()[p.splittingDomain] == minPt) {
      median++;
    }
    if (median == regions.end()) {
      // every possible split in this dimension is the minimum :(
      // so try the next dimension
      p.splittingDomain = (p.splittingDomain + 1) % DIM;
      return calcRegionSplit(p);
    }
    return (*median).box.min()[p.splittingDomain];
  }

  /** Splits the given region page. Will also split parent pages that overflow
      as a result of the split
   * @param[in] p The page to split
   * @param[in] right_parent Optional pointer to page that should be the right page of the split's parent. Set to NULL to make it just p.parent
   * @param[in] dim Dimension to split on 0 <= dim < DIM
   * @param[in] split_pt The coordinate (of @a dim) to split on
   * @pre p.children.size() > n_crit_region_
   * @post p.children.size() <= n_crit_region_
   */
  void split_rp(RegionPage& p,  RegionPage *right_parent, unsigned dim, double split_pt) {
     RegionPage &right_page = *(new RegionPage());
     regionPages_.push_back(&right_page);
     
     p.splittingDomain = right_page.splittingDomain = (p.splittingDomain+1) % DIM;

     // save old (current) regions so we can look through them all
     std::vector<region_type> old_left_regions = p.children;
     // clear the vector so we can add only regions that aren't going into right_page
     p.children = std::vector<region_type>();

     // ensure we have a parent (e.g. if we're currently the root)
     create_parent(&p);
     
     auto cur_box = p.parent->children[p.pidx].box;
     point_type left_max = cur_box.max();
     // left box ends at split_pt in the dim dimension
     left_max[dim] = split_pt;
     point_type right_min = cur_box.min();
     // right box begins at split_pt in the dim dimension
     right_min[dim] = split_pt;
     bounding_box_type left_box(cur_box.min(), left_max);
     bounding_box_type right_box(right_min, cur_box.max());
     
     auto rt1 = region_type{left_box, &p};
     // need to update parent with new box information
     update_parent(rt1);
     if (!right_parent)
       right_parent = p.parent;
     assert(right_parent);

     for (auto&& region : old_left_regions) {
       assert(!(left_box.contains(region.box) && right_box.contains(region.box)));
       if (left_box.contains(region.box)) { // fully left of split
         p.children.push_back(region);
         region.page->pidx = p.children.size()-1;
       } else if (right_box.contains(region.box)) { // fully right of split
         right_page.children.push_back(region);
         region.page->pidx = right_page.children.size()-1;
         region.page->parent = &right_page;
       } else { // middle of split
	 // the current region.page part will stay in our (p's) children
         p.children.push_back(region);
         region.page->pidx = p.children.size()-1;
	 // and the split off right part of region.page goes in right_page's children
         split(*region.page, &right_page, dim, split_pt);
       }
     }
     // we should be all good now
     assert(p.children.size() <= n_crit_region_);
     assert(right_page.children.size() <= n_crit_region_);

     // push the right page we just added to the parent
     auto rt2 = region_type{right_box, &right_page};
     push_page(right_parent, rt2);
   }

  /** Determines where to split given PointPage at
   * (Uses p.splittingDomain as dimension)
   */
  double calcPointSplit(PointPage& p) {
    auto median = p.points.begin() + p.points.size()/2;
    auto dim = p.splittingDomain;
    // median of key_i
    std::nth_element (p.points.begin(), median, p.points.end(),
                      [=] (const point_type& p1, const point_type& p2) {
                        return p1[dim] < p2[dim];
                      });
    return (*median)[p.splittingDomain];
  }
  
  /** Splits the given point page. Will also adjust any parent pages that
      overflow as a result of split
   * @param[in] p The point page to split
   * @param[in] right_parent Optional pointer to page that should be the right page of the split's parent. Set to NULL to make it just p.parent
   * @param[in] dim Dimension to split on 0 <= dim < DIM
   * @param[in] split_pt The coordinate (of @a dim) to split on
   * @pre p.points.size() > n_crit_point_
   * @post p.points.size() <= n_crit_point_
   */
  void split_pp(PointPage &p,  RegionPage *right_parent, unsigned dim, double split_pt) {
    // create new right page, the old page will be left
    PointPage &new_p = *(new PointPage());
    pointPages_.push_back(&new_p);

    // new points for p
    std::vector<point_type> ppoints;

    // init new region_types
    create_parent(&p);
    auto cur_box = p.parent->children[p.pidx].box;
    point_type left_max = cur_box.max();
    // box ends at split_pt in the dim dimension
    left_max[dim] = split_pt;
    point_type right_min = cur_box.min();
    // box begins at split_pt in the dim dimension
    right_min[dim] = split_pt;

    bounding_box_type left_box(cur_box.min(), left_max);
    bounding_box_type right_box(right_min, cur_box.max());

    // reassign the points of the original region
    for (auto&& pnt: p.points) {
      if (pnt[dim] < split_pt)
        ppoints.push_back(pnt);
      else
        new_p.points.push_back(pnt);
    }

    // update p's points and the splittingDomain
    p.points = ppoints;
    p.splittingDomain = (p.splittingDomain + 1) % DIM;
    new_p.splittingDomain = p.splittingDomain;
    
    // update parent to hold new bounding_box
    region_type rt_p = region_type{left_box, &p};    
    update_parent(rt_p);

    if (!right_parent)
      right_parent = p.parent;
    assert(right_parent);

    // add new right_page to its parent
    region_type rt_np = region_type{right_box, &new_p};
    push_page(right_parent, rt_np);
   } 

  /** Update rt.page to have the given region_type in its parent 
   * (Essentially, updating rt.page's box in its parent)
   */
  void update_parent(region_type &rt) {
    assert(rt.page->parent != NULL);
    rt.page->parent->children[rt.page->pidx].box = rt.box;
  }

  /** create a parent for the given page if necessary
   * (i.e. if it is currently the root page)
   */
  void create_parent(Page *page) {
    if (page->parent)
      return;
    RegionPage *parent = new RegionPage();
    regionPages_.push_back(parent);
    parent->children.push_back(region_type{rootBox, page});
    page->parent = parent;
    page->pidx = 0;
    root = parent;
  }
  
  /** Add @a rt into @a parent (and split parent if it then overflows) */
  void push_page(RegionPage *parent, region_type &rt) {
    parent->children.push_back(rt);
    rt.page->parent = parent;
    rt.page->pidx = parent->children.size() - 1;
    if (parent->children.size() > n_crit_region_) {
      split(*parent);
    }
  }
  
  /** Split the given Page
   * @param[in] p Page to split
   * @param[in] right_parent Optional pointer to page that should be the right page of the split's parent. Set to NULL to make it just p.parent (used for pushing the split of a region page downward)
   * @pre p is a full point page the first time this is called
   * @pre left parent is p.parent
   * @post p is in a valid state
   */
  void split(Page &p, RegionPage *right_parent = NULL) {
    double split_pt = p.isRegionPage ? 
    calcRegionSplit((dynamic_cast<RegionPage&> (p))) : calcPointSplit((dynamic_cast<PointPage&> (p)));
    split(p, right_parent, p.splittingDomain, split_pt);
  }

  void split(Page &p, RegionPage *right_parent, unsigned dim, double split_pt) {
      if (p.isRegionPage)
        split_rp(((dynamic_cast<RegionPage&> (p))), right_parent, dim, split_pt);
      else
        split_pp(((dynamic_cast<PointPage&> (p))), right_parent, dim, split_pt);
  }
  
  public:
  /** Returns true iff @a p is in the tree.
   * @param[in] p Point to search for
   * Complexity: O(logn)
   */
  bool query(point_type p) {
    PointPage *page = _query(p);
    return inPointPage(p, page);
  }

private:
  /** Return true iff @a p is in the given PointPage */
  bool inPointPage(point_type p, PointPage *page) {
    for (auto&& point : page->points) {
      if (point == p) {
	return true;
      }
    }
    return false;
  }

  /** Finds the PointPage that a given point would be on
   * @param[in] p Point to search for
   * @param[in] head Pointer to node to start at (usually just used for recursion)
   * @return PointPage that p would be in, if it's in the tree
   */
  PointPage *_query(point_type p, Page *head = NULL) {
    if (!head)
      head = root;
    // we reached the leaf
    if (!head->isRegionPage) {
      return (dynamic_cast<PointPage*> (head));
    }
  
    RegionPage &rp = *(dynamic_cast<RegionPage *> (head));
    for (auto&& region : rp.children) {
      assert(region.page->parent == &rp);
      assert(region.box == region.page->parent->children[region.page->pidx].box);
      // have to use our own contain function here, because we have the added
      // semantics that we only consider a point on the boundary of a box
      // inside that box if it is on the "left" (minimum coordinate) boundary
      bool contained = true;
      auto min = region.box.min();
      auto max = region.box.max();
      for (unsigned i = 0; i != DIM; ++i) {
        if (p[i] < min[i] || p[i] >= max[i]) {
          contained = false;
          break;
        }
      }
      if (contained) {
        return _query(p, region.page);
      }
    }

    assert(0);
    return NULL;
  }

  /** Makes a bounding box containing all of the points in [first, last) */
  template <typename PointIter>
  bounding_box_type get_boundingbox(PointIter first, PointIter last) {
    // Construct a bounding box
    bounding_box_type bb(first, last);
    // Determine the size of the maximum side
    point_type extents = bb.dimensions();
    double max_side = *std::max_element(extents.begin(), extents.end());
    // Make it square and add some wiggle room   TODO: Generalize on square
    return bounding_box_type(bb.center(), (1.0+1e-6) * max_side / 2.0);
  }
  
  Page *root;

  bounding_box_type rootBox;
  
  std::vector<PointPage*> pointPages_;
  std::vector<RegionPage*> regionPages_;
};

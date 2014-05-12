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
#include "MortonCoder.hpp"

#define CACHE_SZ 4096

template <unsigned DIM>
struct NDBTree {

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
  BoundingBox<DIM> box;
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


static constexpr unsigned n_crit_region = CACHE_SZ / (sizeof(region_type)) - CACHE_SZ/100;
static constexpr unsigned n_crit_point = CACHE_SZ / (sizeof(point_type)) - CACHE_SZ/100;

  unsigned n_crit_region_;
  unsigned n_crit_point_;
  // MortonCoder<DIM> coder_;

  template <typename PointIter>
  NDBTree(PointIter first, PointIter last, unsigned ncr_usr = n_crit_region, unsigned ncp_usr = n_crit_point)
      : n_crit_region_(ncr_usr), n_crit_point_(ncp_usr) {
    std::cout << "n_crit_region " << n_crit_region_ << ", n_crit_point " << n_crit_point_ << std::endl;
    
    // create empty root page
    PointPage *pp = new PointPage();
    pointPages.push_back(pp);
    root = pp;
    rootBox = get_boundingbox(first, last);
    
    insert_range(first, last);
  }
  
   void print(Page *p) {
      if (!p->isRegionPage) {
          PointPage *pp = dynamic_cast<PointPage*> (p);

          std::cout << "Point Page with bounding box" << std::endl;
          for (auto&& point: pp->points) {
   //         std::cout << point << std::endl;
          }
   //       std::cout << std::endl;
        return;
      }
      RegionPage *rp = dynamic_cast<RegionPage*> (p);

      for (auto&& region: rp->children) {
        std::cout << "Bounding box: " << region.box << std::endl;
        print(region.page);
      }
  }
  void print() {
    //return;
  //    std::cout << "\n\n\n\n Printing Tree " << std::endl;
//    std::cout << "Root bounding box" << rootBox << std::endl;
    print(root);
//    std::cout << "End of tree\n\n" << std::endl;
  }
  
  private:
  

  //! Uses incremental bucket sorting
  template <typename PointIter>
  void insert_range(PointIter p_first, PointIter p_last) {
    
    for (auto it = p_first; it != p_last; ++it) {
      if (!insert(*it)) {
        // inserted duplicate point // TODO: notify user
        assert(0);
      }
    }
  }
  
  // insert a single point
  bool insert(point_type p) {
    // we create the root in the constructor
    assert(root != NULL);
    assert(rootBox.contains(p));

    //std::cout << "Adding point: " << p << std::endl;
    
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

  void split_rp(RegionPage& p,  RegionPage *right_parent, unsigned dim, double split_pt) {
     RegionPage &right_page = *(new RegionPage());
     
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
     BoundingBox<DIM> left_box(cur_box.min(), left_max);
     BoundingBox<DIM> right_box(right_min, cur_box.max());
     
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

  void split_pp(PointPage &p,  RegionPage *right_parent, unsigned dim, double split_pt) {
    // create new right page, the old page will be left
    PointPage &new_p = *(new PointPage());
    pointPages.push_back(&new_p);

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
    BoundingBox<DIM> left_box(cur_box.min(), left_max);
    BoundingBox<DIM> right_box(right_min, cur_box.max());

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
   
   void update_parent(region_type &rt) {
     assert(rt.page->parent != NULL);
     rt.page->parent->children[rt.page->pidx].box = rt.box;
   }

  void create_parent(Page *page) {
    if (page->parent)
      return;
    RegionPage *parent = new RegionPage();
    parent->children.push_back(region_type{rootBox, page});
    page->parent = parent;
    page->pidx = 0;
    root = parent;
  }
        
   void push_page(RegionPage *parent, region_type &rt) {
    parent->children.push_back(rt);
    rt.page->parent = parent;
    rt.page->pidx = parent->children.size() - 1;
    if (parent->children.size() == n_crit_region_) {
      split(*parent);
    }
   }
  
  /* @pre p is a full point page the first time this is called
   * @pre left parent is p.parent
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
   * @param[out] page will be set to the PointPage that either contains or would contain @a p
   */
  bool query(point_type p) {
    PointPage *page = _query(p);
    return inPointPage(p, page);
  }

private:
  bool inPointPage(point_type p, PointPage *page) {
    for (auto&& point : page->points) {
      if (point == p) {
	return true;
      }
    }
    return false;
  }

  /** finds just the PointPage that a given point would be on */
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

      std::cout << "pidxtwo " << rp.pidx << std::endl;
        std::cout << "My region_box " << rp.parent->children[rp.pidx].box << std::endl;
            std::cout << "My rootbox " << rootBox << std::endl;
    std::cout << "Tried to insert " << p << std::endl;
        assert(rp.parent->children[rp.pidx].box.contains(p));
        if (rp.isRegionPage) {
          for (auto&& pg : rp.children ) {
            std::cout << "Children region box" << pg.box << std::endl;
          }
        }

    std::cout << "My rootbox " << rootBox << std::endl;
    std::cout << "Tried to insert " << p << std::endl;
    // print(root);
    assert(0);
    return false;
  }

  template <typename PointIter>
  BoundingBox<DIM> get_boundingbox(PointIter first, PointIter last) {
    // Construct a bounding box
    BoundingBox<DIM> bb(first, last);
    // Determine the size of the maximum side
    point_type extents = bb.dimensions();
    double max_side = *std::max_element(extents.begin(), extents.end());
    // Make it square and add some wiggle room   TODO: Generalize on square
    return BoundingBox<DIM>(bb.center(), (1.0+1e-6) * max_side / 2.0);
  }
  
  Page *root;

  BoundingBox<DIM> rootBox;
  
  std::vector<PointPage*> pointPages;
};

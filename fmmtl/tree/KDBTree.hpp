#pragma once
/** @file NDTree
 * @brief General class representing a {1D,2D,3D,4D}-Tree.
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


struct Page {
    bool isRegionPage;
    uint8_t splittingDomain;
    uint16_t pidx;
    Page *parent;
};

struct RegionPage : Page {
    struct region_type {
      BoundingBox<DIM> box;
      Page *page;
    };
    
    std::vector<region_type> children;
    RegionPage() : isRegionPage(true), splittingDomain(0), pidx(0), parent(NULL), children() {}

};

struct PointPage : Page {
    std::vector<point_type> points;
    PointPage() : isRegionPage(false), splittingDomain(0), pidx(0), parent(NULL), points() {}
};


static constexpr unsigned n_crit_region = CACHE_SZ / (sizeof(RegionPage::region_type));
static constexpr unsigned n_crit_point = CACHE_SZ / (sizeof(point_type));

  unsigned n_crit_region_;
  unsigned n_crit_point_;
  MortonCoder<DIM> coder_;

  template <typename PointIter>
  NDBTree(PointIter first, PointIter last, unsigned ncr_usr = n_crit_region, unsigned ncp_usr = n_crit_point)
      : coder_(get_boundingbox(first, last)), n_crit_region_(ncr_usr), n_crit_point_(ncp_usr) {

    // create empty root page
    pointPages.emplace_back();
    root = &pointPages[0];
    
    insert_range(first, last, n_crit_point_usr, n_crit_point_region_usr);
  }
  
   void print() {
      for (auto&& p: pointPages) {
        std::cout << "Point Page \n";
        for (auto&& point: p.points) {
          std::cout << point << ", ";
        }
        std::cout << std::endl;
      }
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
    
    PointPage page;
    if (query(p, page)) {
      return false;
    }
    
    page.points.push_back(p);
    
    if (page.points.size() <= n_crit_point_) {
      return true;
    }
    
    split(page);
    
    return true;
  }

  double calcRegionSplit(const RegionPage& p) {
    auto& regions = p.children;
    auto median = regions.begin() + regions.size()/2;
    std::nth_element(regions.begin(), median, regions.end(), 
                     [] (region_type r1, region_type r2) { return r1.box.min()[p.splittingDomain] < r2.box.min()[p.splittingDomain]; });
    return (*median).box.min()[p.splittingDomain];
  }

   void split_rp(RegionPage& p,  Page *right_parent = p.parent) {     
     regionPages.emplace_back();
     RegionPage& right_page = regionPages.back();
     
     p.splittingDomain = right_page.splittingDomain = (p.splittingDomain + 1) % DIM;

     // save old (current) regions so we can look through them all
     std::vector<region_type> old_left_regions = p.children;
     // clear the vector so we can add only regions that aren't going into right_page
     p.children = std::vector<region_type>();
     
     auto cur_box = p.parent->children[p.pidx].box;
     point_type left_max = cur_box.max();
     // box ends at split_pt in the dim dimension
     left_max[dim] = split_pt;
     point_type right_min = cur_box.min();
     // box begins at split_pt in the dim dimension
     right_min[dim] = split_pt;
     BoundingBox<DIM> left_box(cur_box.min(), left_max);
     BoundingBox<DIM> right_box(right_min, cur_box.max());
     
     update_parent(region_type{left_box, &p});
     
     for (auto&& region : old_left_regions) {
       if (left_box.contains(region.box)/*region.box.max()[dim] <= split_pt*/) { // fully left of split
         p.children.push_back(region);
       } else if (right_box.contains(region.box)/*region.box.min()[dim] >= split_pt*/) { // fully right of split
         right_page.children.push_back(region);
       } else { // middle of split
         //TODO: refactor split to take dim and splitpt
         split(region.page, &right_page, dim, split_pt);
       }
     }
     p.parent->children[p.pidx].box = left_box;
     push_page(region_type{right_box, &right_page}, p->parent);
   }

   double calcPointSplit(const PointPage& p) {
    auto median = p.points.begin() + p.points.size()/2;
    // median of key_i
    std::nth_element (p.points.begin(), median, p.points.end(),
                      [] (point_type p1, point_type p2) {
                        return p1[p.splittingDomain] < p2[p.splittingDomain];
                      });
    return (*median)[p.splittingDomain];
  }

   void split_pp(PointPage p,  Page *right_parent = p.parent) {
    // median of key_i
    std::nth_element (p.points.begin(), p.points.end(), 
      [] (point_type p1, point_type p2) { 
        return p1[p.splittingDomain] < p2[p.splittingDomain];
      });
    point_type median = p.points[p.points.size() / 2];

    // create new right page, the old page will be left
    pointPages_.emplace_back();
    PointPage &new_p = pointPages.back();
    // new points for p
    std::vector<point_type> ppoints;

    // init new region_types
    RegionPage::region_type rt_p = {.box = BoundingBox(), .page = p};
    RegionPage::region_type rt_np = {.box = BoundingBox(), .page = &new_page};

    // reassign the points of the original region
    for (auto&& pnt: p.points) {
      if (pnt[p.splittingDomain] < median) {
        ppoints.push_back(pnt);
        rt_p.box |= pnt; //
      }
      else {
        new_p.points.push_back(pnt);
        rt_np.box |= pnt;
      }
    }
    // update p's points and the splittingDomain
    p.points = ppoints;
    p.splittingDomain = (p.splittingDomain + 1) % n_crit_point_;
    new_p.splittingDomain = p.splittingDomain;
    
    // update parent to hold new bounding_box
    update_parent(rt_p);
    return push_page(p->parent, rt_np);
   } 
   
   void update_parent(region_type &rt) {
      if (rt.p->parent != NULL) {
         p.parent->children[p.pid].box = rt.box;
      }
      else { // a parent does not exist so create one
        regionPages_.emplace_back();
        RegionPage &parent = regionPages.back();
        parent.children.push_back(rt);
        rt.p->parent = parent;
        rt.p->pid = 0;
        root = &parent;
      }
   }
        
   void push_page(RegionPage *parent, region_type &rt) {
    parent->children.push_back(rt);
    rt.page->parent = parent;
    rt.page->pid = parent->children.size() - 1;
    if (parent->children.size() == n_crit_region_) {
      split(*parent);
    }
   }
  
  /* @pre p is a full point page the first time this is called
   * @pre left parent is p.parent
   */

  void split(Page &p, Page *right_parent = NULL) {
    if (!right_parent) {
      right_parent = p.parent;
    }
    double split_pt = p.isRegionPage ? calcRegionSplit(p) : calcPointSplit(p);
    split(p, right_parent, p.splittingDomain, split_pt);
  }

  void split(Page &p, Page *right_parent, unsigned dim, double split_pt) {
      if (p.isRegionPage)
        split_rp(p, right_parent, dim, split_pt);
      else
        split_pp(p, right_parent, dim, split_pt);
  }
  
  /** Returns true iff @a p is in the tree.
   * @param[out] page will be set to the PointPage that either contains or would contain @a p
   */
  bool query(point_type p, PointPage& pointpage, Page *head = root) {
    // we reached the leaf
    if (!head->isRegionPage) {
      pointpage = head;
      // TODO: we could speed this up slightly if we just ignored/assumed no duplicate points
      // (by not doing this for loop)
      for (auto&& point : pointpage.points) {
          if (point == p) {
              return true;
          }
      }
      return false;
    }
  
    for (auto&& region : head->children) {
      if (region.box.contains(p)) {
        return query(p, pointpage, region.page);
      }
    }
  }
  

  Page *root;
  
  std::vector<PointPage> pointPages;
  std::vector<RegionPage> regionPages;
};

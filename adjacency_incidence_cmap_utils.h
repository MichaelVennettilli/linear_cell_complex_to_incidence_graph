/*
MIT License

Copyright (c) 2022 Michael Vennettilli

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

#include <vector>
#include <random>
#include <math.h>
#include <cmath>
#include <map>
#include <utility>

// Typedefs pertaining to the convex hull LCC
typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_CH;
typedef LCC_CH::typename CMap::Dart_handle           typename CMap::Dart_handle_CH;
// Typedefs pertaining to the shell LCC
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
typedef LCC_3::typename CMap::Dart_handle           typename CMap::Dart_handle_3;
typedef LCC_3::Vertex_attribute_handle           Vertex_handle_3;

// Inverse function for indices
inline unsigned int inv(unsigned int i){
  return (i>1) ? i : (1-i);
}

// Mark all darts in c_i(dh)
template<typename CMap>
void mark_darts_in_cell(CMap &cm, unsigned int i, typename CMap::Dart_handle &dh, typename CMap::size_type ma){
  for(typename CMap::Dart_of_cell_range<i>::iterator
    dart_it=cm.darts_of_cell<i>(dh).begin(),
    dart_it_end=cm.darts_of_cell<i>(dh).end();
    dart_it!=dart_it_end; ++dart_it){
    cm.mark(dart_it, ma);
  }
  return;
}

// Unmark all darts in c_i(dh)
template<typename CMap>
void unmark_darts_in_cell(CMap &cm, unsigned int i, typename CMap::Dart_handle &dh, typename CMap::size_type ma){
  for(typename CMap::Dart_of_cell_range<i>::iterator
    dart_it=cm.darts_of_cell<i>(dh).begin(),
    dart_it_end=cm.darts_of_cell<i>(dh).end();
    dart_it!=dart_it_end; ++dart_it){
    cm.unmark(dart_it, ma);
  }
  return;
}

// Check if a dart in c_i(dh) has already been marked
template<typename CMap>
bool is_dart_in_cell_marked(CMap &cm, unsigned int i, typename CMap::Dart_handle &dh, typename CMap::size_type ma){
  bool result = false;
  typename CMap::Dart_of_cell_range<i>::iterator dart_it=cm.darts_of_cell<i>(dh).begin();
  typename CMap::Dart_of_cell_range<i>::iterator dart_it_end=cm.darts_of_cell<i>(dh).end();
  while((!result)&&(dart_it!=dart_it_end)){
    if(cm.is_marked(dart_it,ma)) result = true;
    ++dart_it;
  }
  return result;
}

// Returns true if the i-cell corresponding to dh_i is incident to the
// j-cell corresponding to dh_j.
template<typename CMap>
bool check_incident_pair(CMap &cm, unsigned int i, typename CMap::Dart_handle &dh_i,
  unsigned int j, typename CMap::Dart_handle &dh_j){
  typename CMap::size_type ma = cm.get_new_mark();
  mark_darts_in_cell(cm, i, dh_i, ma);
  bool result = is_dart_in_cell_marked(cm, j, dh_j, ma);
  unmark_darts_in_cell(cm, i, dh_i, ma);
  cm.free_mark(ma);
  return result;
}

// Returns true if c_i(dh_1) is adjacent to c_i(dh_2). Work with the case i > 0 first.
template<typename CMap>
bool check_adjacent_pair(CMap &cm, unsigned int i, typename CMap::Dart_handle &dh_1, typename CMap::Dart_handle &dh_2){
  bool result = false;
  typename CMap::size_type ma = cm.get_new_mark();
  mark_darts_in_cell(cm, i, dh_1, ma);

  typename CMap::Dart_of_cell_range<i>::iterator dart_it=cm.darts_of_cell<i>(dh_2).begin();
  typename CMap::Dart_of_cell_range<i>::iterator dart_it_end=cm.darts_of_cell<i>(dh_2).end();
  if(i==0){
    unsigned int max_dimension = cm.dimension;
    while((!result)&&(dart_it!=dart_it_end)){
      for(unsigned int dim_it = 0; dim_it <= max_dimension; dim_it++){
        if(cm.is_marked(cm.beta(dart_it,i),ma)) result = true;
      }
      ++dart_it;
    }
  } else {
    // Loop over things in c_i(dh_2) and see if applying Beta_i  of Beta_{inv(i)}
    // gets you to something that has already been hit.
    while((!result)&&(dart_it!=dart_it_end)){
      if(cm.is_marked(cm.beta(dart_it,i),ma)) result = true;
      if(cm.is_marked(cm.beta(dart_it, inv(i)),ma)) result = true;
      ++dart_it;
    }
  }
  // Unmark everything and free the mark
  unmark_darts_in_cell(cm, i, dh_1, ma);

  cm.free_mark(ma);
  return result;
}

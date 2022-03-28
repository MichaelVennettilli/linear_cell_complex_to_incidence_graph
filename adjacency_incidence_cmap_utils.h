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
unsigned int inv(unsigned int i){
  return (i>1) ? i : (1-i);
}

// Returns true if the i-cell corresponding to dh_i is incident to the
// j-cell corresponding to dh_j.
template<typename CMap, unsigned int i, unsigned int j>
bool check_incident_pair(CMap &cm, typename CMap::Dart_handle &dh_i, typename CMap::Dart_handle &dh_j){
  bool result = false;
  typename CMap::size_type ma = cm.get_new_mark();
  // Mark everything in the i-cell belonging to c_i(dh_i)
  for(typename CMap::Dart_of_cell_range<i>::iterator
    i_cell_it=cm.darts_of_cell<i>(dh_i).begin(),
    i_cell_end=cm.darts_of_cell<i>(dh_i).end();
    i_cell_it!=i_cell_end; ++i_cell_it){
    cm.mark(i_cell_it, ma);
  }
  // Loop over darts in c_j(dh_j) until you reach the end or hit something marked
  typename CMap::Dart_of_cell_range<j>::iterator j_cell_it=cm.darts_of_cell<j>(dh_j).begin();
  typename CMap::Dart_of_cell_range<j>::iterator j_cell_end=cm.darts_of_cell<j>(dh_j).end();
  while((!result)&&(j_cell_it!=j_cell_end)){
    if(cm.is_marked(j_cell_it,ma)) result = true;
    ++j_cell_it;
  }
  // Unmark everything and free the mark
  for(typename CMap::Dart_of_cell_range<i>::iterator
    i_cell_it=cm.darts_of_cell<i>(dh_i).begin(),
    i_cell_end=cm.darts_of_cell<i>(dh_i).end();
    i_cell_it!=i_cell_end; ++i_cell_it){
    cm.unmark(i_cell_it, ma);
  }
  cm.free_mark(ma);
  return result;
}

// Returns true if c_i(dh_1) is adjacent to c_i(dh_2). Work with the case i > 0 first.
template<typename CMap, unsigned int i>
bool check_adjacent_pair(CMap &cm, typename CMap::Dart_handle &dh_1, typename CMap::Dart_handle &dh_2){
  bool result = false;
  typename CMap::size_type ma = cm.get_new_mark();
  // Mark everything in c_i(dh_1)
  for(typename CMap::Dart_of_cell_range<i>::iterator
    i_cell_it=cm.darts_of_cell<i>(dh_1).begin(),
    i_cell_end=cm.darts_of_cell<i>(dh_1).end();
    i_cell_it!=i_cell_end; ++i_cell_it){
    cm.mark(i_cell_it, ma);
  }
  // Loop over things in c_i(dh_2) and see if applying Beta_i gets you
  // to something that has already been hit.
  typename CMap::Dart_of_cell_range<i>::iterator i_cell_it=cm.darts_of_cell<i>(dh_2).begin();
  typename CMap::Dart_of_cell_range<i>::iterator i_cell_end=cm.darts_of_cell<i>(dh_2).end();
  while((!result)&&(i_cell_it!=i_cell_end)){
    if(cm.is_marked(cm.beta(i_cell_it,i),ma)) result = true;
    if(cm.is_marked(cm.beta(i_cell_it, inv(i)),ma)) result = true;
    ++j_cell_it;
  }
  // Unmark everything and free the mark
  for(typename CMap::Dart_of_cell_range<i>::iterator
    i_cell_it=cm.darts_of_cell<i>(dh_1).begin(),
    i_cell_end=cm.darts_of_cell<i>(dh_1).end();
    i_cell_it!=i_cell_end; ++i_cell_it){
    cm.unmark(i_cell_it, ma);
  }
  cm.free_mark(ma);
  return result;
}

// Returns true if c_0(dh_1) is adjacent to c_0(dh_2)
template<typename CMap, unsigned int i>
bool check_adjacent_pair<CMap, 0>(CMap &cm, typename CMap::Dart_handle &dh_1, typename CMap::Dart_handle &dh_2){
  bool result = false;
  typename CMap::size_type ma = cm.get_new_mark();
  int dim = cm.dimension;
  // Mark everything in c_i(dh_1)
  for(typename CMap::Dart_of_cell_range<i>::iterator
    dart_it=cm.darts_of_cell<0>(dh_1).begin(),
    dart_it_end=cm.darts_of_cell<0>(dh_1).end();
    dart_it!=dart_it_end; ++dart_it){
    cm.mark(dart_it, ma);
  }
  // Loop over things in c_i(dh_2) and see if applying Beta_i gets you
  // to something that has already been hit.
  typename CMap::Dart_of_cell_range<i>::iterator dart_it=cm.darts_of_cell<0>(dh_2).begin();
  typename CMap::Dart_of_cell_range<i>::iterator dart_it_end=cm.darts_of_cell<0>(dh_2).end();
  while((!result)&&(dart_it!=dart_it_end)){
    for(int i=0; i<= dim; i++){
      if(cm.is_marked(cm.beta(dart_it,i),ma)) result = true;
    }
    ++dart_it;
  }
  // Unmark everything and free the mark
  for(typename CMap::Dart_of_cell_range<i>::iterator
    dart_it=cm.darts_of_cell<0>(dh_1).begin(),
    dart_it_end=cm.darts_of_cell<0>(dh_1).end();
    dart_it!=dart_it_end; ++dart_it){
    cm.unmark(dart_it, ma);
  }
  cm.free_mark(ma);
  return result;
}

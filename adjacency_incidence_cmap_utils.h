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

#ifndef ADJACENCY_INCIDENCE_CMAP_UTILS_H
#define ADJACENCY_INCIDENCE_CMAP_UTILS_H

// Compute the index of the inverse beta transformation.
unsigned int inv(unsigned int i){
  return (i>1) ? i : (1-i);
}

// Mark all darts in a cell
template<typename CMap, unsigned int cell_dim>
void mark_darts_in_cell(CMap &cm, typename CMap::Dart_handle &dh, typename CMap::size_type ma){
  auto darts_of_cell_range = cm.template darts_of_cell<cell_dim>(dh);
  for(auto dart_it = darts_of_cell_range.begin(),
  dart_it_end = darts_of_cell_range.end(); dart_it!=dart_it_end; dart_it++){
    cm.mark(dart_it, ma);
  }
  return;
}

// Unmark all darts in a cell
template<typename CMap, unsigned int cell_dim>
void unmark_darts_in_cell(CMap &cm, typename CMap::Dart_handle &dh, typename CMap::size_type ma){
  auto darts_of_cell_range = cm.template darts_of_cell<cell_dim>(dh);
  for(auto dart_it = darts_of_cell_range.begin(),
  dart_it_end = darts_of_cell_range.end(); dart_it!=dart_it_end; dart_it++){
    cm.unmark(dart_it, ma);
  }
  return;
}

// Check if some dart in a cell has been marked
template<typename CMap, unsigned int cell_dim>
bool is_dart_in_cell_marked(CMap &cm, typename CMap::Dart_handle &dh, typename CMap::size_type ma){
  bool result = false;
  auto darts_of_cell_range = cm.template darts_of_cell<cell_dim>(dh);
  auto dart_it = darts_of_cell_range.begin();
  auto dart_it_end = darts_of_cell_range.end();
  while((!result)&&(dart_it!=dart_it_end)){
    if(cm.is_marked(dart_it,ma)) result = true;
    dart_it++;
  }
  return result;
}

// Check if the cells c_i(dh_i) and c_j(dh_j) are incident
template<typename CMap, unsigned int i, unsigned int j>
bool check_incident_pair(CMap &cm, typename CMap::Dart_handle &dh_i, typename CMap::Dart_handle &dh_j){
  typename CMap::size_type ma = cm.get_new_mark();
  mark_darts_in_cell<CMap,i>(cm,dh_i,ma);
  bool result = is_dart_in_cell_marked<CMap,j>(cm,dh_j,ma);
  unmark_darts_in_cell<CMap,i>(cm,dh_i,ma);
  cm.free_mark(ma);
  return result;
}

// Check if the cells corresponding to dh_1 and dh_2 are adjacent.
template<typename CMap, unsigned int cell_dim>
bool check_adjacent_pair(CMap &cm, typename CMap::Dart_handle &dh_1, typename CMap::Dart_handle &dh_2){
  bool result;
  typename CMap::size_type ma = cm.get_new_mark();
  auto darts_of_cell_range = cm.template darts_of_cell<cell_dim>(dh_2);
  auto dart_it = darts_of_cell_range.begin();
  auto dart_it_end = darts_of_cell_range.end();

  mark_darts_in_cell<CMap,cell_dim>(cm,dh_1,ma);
  if(cell_dim>0){
    while((!result)&&(dart_it!=dart_it_end)){
      if(cm.is_marked(cm.beta(dart_it,cell_dim),ma)) result = true;
      if(cm.is_marked(cm.beta(dart_it,inv(cell_dim)),ma)) result = true;
      ++dart_it;
    }
  } else {
    unsigned int cmap_dim = cm.dimension;
    while((!result)&&(dart_it!=dart_it_end)){
      for(int i=0; i<= cmap_dim; i++){
        if(cm.is_marked(cm.beta(dart_it,i),ma)) result = true;
      }
    }
  }
  unmark_darts_in_cell<CMap,cell_dim>(cm,dh_1,ma);
  cm.free_mark(ma);
  return result;
}

#endif

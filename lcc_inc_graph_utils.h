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

#include "adjacency_incidence_cmap_utils.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <utility>


typedef std::vector<std::pair<std::string, std::string>> String_edge_list;

template<typename Lcc>
void write_inc(Lcc &lcc){
  String_edge_list edge_list;
  std::vector<std::map<typename Lcc::Dart_handle,std::string>> dictionary_vector[lcc.dimension+1];
  int counter = 0;
  typename Lcc::size_type ma = lcc.get_new_mark();

  // Add the entries to the list of dictionaries
  for(int i=0; i<=lcc.dimension; i++){
    counter = 0;
    for(typename Lcc::template One_dart_per_cell_range<i>::iterator
      dart_it=lcc.one_dart_per_cell<i>().begin(), dart_it_end=lcc.one_dart_per_cell<i>().end();
      dart_it!=dart_it_end;++dart_it){
      dictionary_vector[i][dart_it] = std::to_string(i) + "_" + std::to_string(counter);
    }
  }
  // Determine incidence and write the pairs
  for(int i=lcc.dimension; i>0;i--){
    // Iterate over the high dimensional cells
    for(typename Lcc::template One_dart_per_cell_range<i>::iterator
      high_it=lcc.one_dart_per_cell<i>().begin(), high_end=lcc.one_dart_per_cell<i>().end();
      high_it!=high_end;++high_it){
      // Mark the given cell, then loop over lower dimensional cells.
      mark_darts_in_cell<Lcc,i>(lcc, high_it, ma);
      for(typename Lcc::template One_dart_per_cell_range<i-1>::iterator
        low_it=lcc.one_dart_per_cell<i-1>().begin(), low_end=lcc.one_dart_per_cell<i-1>().end();
        low_it!=low_end;++low_end){
        if(check_incident_pair<Lcc,i,i-1>(Lcc, high_it, low_it)){
          edge_list.push_back(make_pair(dictionary_vector[i][high_it], dictionary_vector[i-1][low_it]));
        }
      }
      unmark_darts_in_cell<Lcc,i>(lcc, high_it, ma);
    }
  }
  for(int i=0; i< edge_list.size();i++){
    std::cout << edge_list[i].first << ", " << edge_list[i].second << std::endl;
  }
}

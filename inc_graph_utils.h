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

#ifndef INC_GRAPH_UTILS_H
#define INC_GRAPH_UTILS_H

#include "adjacency_incidence_cmap_utils.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <fstream>


/*
 The incidence graph is stored as an edge list. Here, the geometric elements are
 viewed as nodes in a graph. Nodes are connected by an edge if they are incident/.
*/
typedef std::vector<std::pair<std::string, std::string>> Incidence_graph_by_layers;

template<typename CMap, int max_dimension>
void export_incidence_graph(CMap &cm, std::string filename){
  // Create aliases for the integer sequences
  using seq_for_maps = std::make_integer_sequence<int, max_dimension + 1>;
  using seq_for_compare = std::make_integer_sequence<int, max_dimension>;
  // Create the vector of maps for converting darts to strings.
  std::vector<std::map<typename CMap::Dart_handle,std::string>> dictionary_vector;
  // Create the incidence graph.
  Incidence_graph_by_layers incidence_graph;
  // Reserve a mark
  typename CMap::size_type ma = cm.get_new_mark();

  // Create the lambda that will sets the i-th component of the dictionary vector.
  auto make_dictionary_layer = [&]<int i>(std::integral_constant<int, i> dim){
    int counter = 0;
    std::map<typename CMap::Dart_handle,std::string> dictionary_layer;
    for(typename CMap::template One_dart_per_cell_range<dim>::iterator
      dart_it=cm.template one_dart_per_cell<dim>().begin(),
      dart_it_end=cm.template one_dart_per_cell<dim>().end();
      dart_it!=dart_it_end;++dart_it){
      dictionary_layer[dart_it] = std::to_string(dim) + "_" + std::to_string(counter++);
    }
    dictionary_vector.push_back(dictionary_layer);
  };

  // Create the lambda that adds one layer to the incidence graph
  auto make_incidence_layer = [&]<int i>(std::integral_constant<int, i> lower_dim,
    std::integral_constant<int, i+1> higher_dim){
      for(typename CMap::template One_dart_per_cell_range<higher_dim>::iterator
        high_it=cm.template one_dart_per_cell<higher_dim>().begin(), high_end=cm.template one_dart_per_cell<higher_dim>().end();
        high_it!=high_end;++high_it){
        // Mark the given cell, then loop over lower dimensional cells.
        std::cout << "marking" << std::endl;
        mark_darts_in_cell<CMap,higher_dim>(cm, high_it, ma);
        std::cout << "marked" << std::endl;
        for(typename CMap::template One_dart_per_cell_range<lower_dim>::iterator
          low_it=cm.template one_dart_per_cell<lower_dim>().begin(),
          low_end=cm.template one_dart_per_cell<lower_dim>().end();
          low_it!=low_end;++low_end){
          std::cout << "check incidence" << std::endl;
          if(check_incident_pair<CMap,higher_dim, lower_dim>(cm, high_it, low_it)){
            incidence_graph.push_back(std::make_pair(dictionary_vector[higher_dim][high_it], dictionary_vector[lower_dim][low_it]));
            std::cout << dictionary_vector[higher_dim][high_it] << ", " << dictionary_vector[lower_dim][low_it] << std::endl;
          }
        }
        std::cout << "unmarking" << std::endl;
        unmark_darts_in_cell<CMap,higher_dim>(cm, high_it, ma);
      }
  };

  // Execute the lambda that iterates over all dictionary layers.
  [&]<int... Is>(std::integer_sequence<int, Is...>){
    (make_dictionary_layer(std::integral_constant<int, Is>{}), ...);
  }(seq_for_maps{});

  // Create the lambda that generates the incidence graph.
  [&]<int... Is>(std::integer_sequence<int, Is...>){
    (make_incidence_layer(std::integral_constant<int, Is>{}, std::integral_constant<int, Is+1>{}), ...);
  }(seq_for_compare{});

  // Write to file
  std::ofstream strm;
  strm.open(filename,std::ofstream::out|std::ofstream::trunc);
  for(int edge_index = 0; edge_index < incidence_graph.size(); edge_index++){
    strm << "(" << incidence_graph[edge_index].first << ", " << incidence_graph[edge_index].second << std::endl;
  }
  return;
}

#endif

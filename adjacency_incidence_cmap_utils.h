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

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <tuple>

#include <fstream>

// Compute the index of the inverse beta transformation.
inline unsigned int inv(unsigned int i){
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

// Check if the i-cells corresponding to dh_1 and dh_2 are adjacent.
template<typename CMap, unsigned int cell_dim>
bool check_adjacent_pair(CMap &cm, typename CMap::Dart_handle &dh_1, typename CMap::Dart_handle &dh_2){
  bool result = false;
  typename CMap::size_type ma = cm.get_new_mark();
  auto darts_of_cell_range = cm.template darts_of_cell<cell_dim>(dh_2);
  auto dart_it = darts_of_cell_range.begin();
  auto dart_it_end = darts_of_cell_range.end();

  mark_darts_in_cell<CMap,cell_dim>(cm,dh_1,ma);
  if(cell_dim>1){
	// In this case, cell_dim == inv(cell_dim)
	while((!result)&&(dart_it!=dart_it_end)){
      if(cm.is_marked(cm.beta(dart_it,cell_dim),ma)) result = true;
      ++dart_it;
    } 
  } else if(cell_dim==1){
    while((!result)&&(dart_it!=dart_it_end)){
      if(cm.is_marked(cm.beta(dart_it,1),ma)) result = true;
      if(cm.is_marked(cm.beta(dart_it,0),ma)) result = true;
      ++dart_it;
    }
  } else {
	// The rule is different for cell_dim == 0
    unsigned int cmap_dim = cm.dimension;
    while((!result)&&(dart_it!=dart_it_end)){
      for(int i=0; i<= cmap_dim; i++){
        if(cm.is_marked(cm.beta(dart_it,i),ma)) result = true;
      }
      ++dart_it;
    }
  }
  unmark_darts_in_cell<CMap,cell_dim>(cm,dh_1,ma);
  cm.free_mark(ma);
  return result;
}

/*
 * The incidence graph is stored as an edge list. Here, the geometric elements are
 * viewed as nodes in a graph. Nodes are connected by an edge if they are incident/.
 */
 
// Helper functions for converting an integer sequence to a parameter
// pack via type deduction.
template<typename CMap, int... Is>
auto make_dict_tuple_helper( std::integer_sequence<int,Is...> ){
    return std::tuple<std::map<typename CMap::template Attribute_handle<std::integral_constant<int, Is>{}>::type, std::string> ...>{};
} 

template<typename CMap, int max>
using Dict_tuple = decltype(make_dict_tuple_helper<CMap>(std::make_integer_sequence<int, max>{}));

template<typename CMap, int max_dimension>
void export_incidence_graph(CMap &cm, std::string filename){
	// Create aliases for the integer sequences
	using seq_for_maps = std::make_integer_sequence<int, max_dimension + 1>;
	using seq_for_compare = std::make_integer_sequence<int, max_dimension>;
	
	// Create the tuple of maps for converting attribute handles to strings.
	// You don't know if getting one dart per cell vs one dart per incident
	// cell will give you the same dart, so you can't use dart handles.
	Dict_tuple<CMap, max_dimension + 1> dict_tuple{};
	// Create the incidence graph.
	typedef std::vector<std::pair<std::string, std::string>> Incidence_graph_by_layers;
	Incidence_graph_by_layers incidence_graph;
	// Reserve a mark
	typename CMap::size_type ma = cm.get_new_mark();
	// Create the printing objects
	std::ofstream strm;
	strm.open(filename,std::ofstream::out|std::ofstream::trunc);
	// Create the lambda that will sets the i-th component of the dictionary vector.
	auto make_dictionary_layer = [&]<int i>(std::integral_constant<int, i> dim){
		std::string as;
		typename CMap::template Attribute_handle<dim>::type ah;
		// Start a counter for the cell number, a proxy for id.
		int counter = 0;
		// Associate each cell to a string
		std::map<typename CMap::Dart_handle,std::string> dictionary_layer;
		for(typename CMap::template One_dart_per_cell_range<dim>::iterator
			dart_it=cm.template one_dart_per_cell<dim>().begin(),
			dart_it_end=cm.template one_dart_per_cell<dim>().end();
			dart_it!=dart_it_end;++dart_it){
			as = std::to_string(dim) + "_" + std::to_string(counter++);
			ah = cm.template attribute<dim>(dart_it);
			std::get<dim>(dict_tuple).insert({ah, as});
		}
		// The first line of the file has the number of cells of each dimension.
		strm << counter;
		if(dim < max_dimension){
			strm << ", ";
		} else {
			strm << std::endl;
		}
	};
	
	// Create the lambda that adds one layer to the incidence graph
	auto make_incidence_layer = [&]<int i>(std::integral_constant<int, i> lower_dim,
		std::integral_constant<int, i+1> higher_dim){
		// Declare strings and create the range of all higher-dimensional cells.
		std::string higher_string, lower_string;
		auto higher_dim_range = cm.template one_dart_per_cell<higher_dim>();
		for(auto high_it=higher_dim_range.begin(), high_end=higher_dim_range.end();
		high_it != high_end; high_it++){
			// Set the string for the higher-D cell and the range of incident lower-D cells
			higher_string = std::get<higher_dim>(dict_tuple)[cm.template attribute<higher_dim>(high_it)];
			auto lower_dim_range = cm.template one_dart_per_incident_cell<lower_dim, higher_dim>(high_it);
			for(auto low_it = lower_dim_range.begin(), low_end = lower_dim_range.end();
			low_it != low_end; low_it++){
				// Get the lower-D string and push the pair to the incidence graph
				lower_string = std::get<lower_dim>(dict_tuple)[cm.template attribute<lower_dim>(low_it)];
				incidence_graph.push_back(std::make_pair(higher_string, lower_string));
			}
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

	// Write the dictionary to file
	for(int edge_index = 0; edge_index < incidence_graph.size()-1; edge_index++){
		strm << incidence_graph[edge_index].first << ", " << incidence_graph[edge_index].second << std::endl;
	}
	strm << incidence_graph[incidence_graph.size()-1].first << ", " << incidence_graph[incidence_graph.size()-1].second;
	cm.free_mark(ma);

	return;
}

// Function for exporting the points of a linear cell complex
template<typename Lcc>
	void export_vertices(Lcc &lcc, std::string filename){
	int dimension = lcc.ambient_dimension;
	// Create storage for the point and counter
	typename Lcc::Point current_point;
	int counter =0;
	// Create the file stream
	std::ofstream strm;
	strm.open(filename,std::ofstream::out|std::ofstream::trunc);
	for(typename Lcc::template One_dart_per_cell_range<0>::iterator
	dart_it=lcc.template one_dart_per_cell<0>().begin(),
	dart_it_end=lcc.template one_dart_per_cell<0>().end();
	dart_it!=dart_it_end;){
		current_point = lcc.point(dart_it);
		for(int coordinate=0; coordinate < dimension-1;coordinate++) strm << current_point[coordinate] << ", ";
		strm << current_point[dimension-1];
		if(++dart_it!=dart_it_end) strm << std::endl;
	}
}

#endif

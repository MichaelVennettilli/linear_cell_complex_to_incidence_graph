template<typename Container, int max_dimension>
void my_function(Container &c){
  using seq = std::make_integer_sequence<int, max_dimension>;
  auto fn_body = [&]<int i>(std::integral_constant<int, i>){
    std::cout << i << std::endl;
  };
  [&]<int... Is>(std::integer_sequence<int, Is...>){
      (fn_body(std::integral_constant<int, Is>{}),...);
  }(seq{});
}

// Old approach

template<typename Lcc, int dimension>
void write_inc(Lcc &lcc){
  Incidence_graph_by_layers edge_list;
  std::vector<std::map<typename Lcc::Dart_handle,std::string>> dictionary_vector[dimension+1];
  int counter = 0;
  typename Lcc::size_type ma = lcc.get_new_mark();

  // Add the entries to the list of dictionaries
  for(int i=0; i<=dimension; i++){
    counter = 0;
    for(typename Lcc::template One_dart_per_cell_range<i>::iterator
      dart_it=lcc.template one_dart_per_cell<i>().begin(), dart_it_end=lcc.template one_dart_per_cell<i>().end();
      dart_it!=dart_it_end;++dart_it){
      dictionary_vector[i][dart_it] = std::to_string(i) + "_" + std::to_string(counter);
    }
  }
  // Determine incidence and write the pairs
  for(int i=dimension; i>0;i--){
    // Iterate over the high dimensional cells
    for(typename Lcc::template One_dart_per_cell_range<i>::iterator
      high_it=lcc.template one_dart_per_cell<i>().begin(), high_end=lcc.template one_dart_per_cell<i>().end();
      high_it!=high_end;++high_it){
      // Mark the given cell, then loop over lower dimensional cells.
      mark_darts_in_cell<Lcc,i>(lcc, high_it, ma);
      for(typename Lcc::template One_dart_per_cell_range<i-1>::iterator
        low_it=lcc.template one_dart_per_cell<i-1>().begin(), low_end=lcc.template one_dart_per_cell<i-1>().end();
        low_it!=low_end;++low_end){
        if(check_incident_pair<Lcc,i,i-1>(Lcc, high_it, low_it)){
          edge_list.push_back(std::make_pair(dictionary_vector[i][high_it], dictionary_vector[i-1][low_it]));
        }
      }
      unmark_darts_in_cell<Lcc,i>(lcc, high_it, ma);
    }
  }
  for(int i=0; i< edge_list.size();i++){
    std::cout << edge_list[i].first << ", " << edge_list[i].second << std::endl;
  }
}

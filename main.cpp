#include <CGAL/Combinatorial_map.h>
#include "adjacency_incidence_cmap_utils.h"
#include "lcc_inc_graph_utils.h"

#include <iostream>

typedef CGAL::Combinatorial_map<2> CM_2;
typedef CM_2::Dart_handle Dart_handle;

int main(){
  CM_2 cm;
  Dart_handle dh1 = cm.create_dart();
  Dart_handle dh2 = cm.create_dart();
  Dart_handle dh3 = cm.create_dart();
  Dart_handle dh4 = cm.create_dart();
  cm.sew<1>(dh1,dh2);
  cm.sew<1>(dh2,dh3);
  cm.sew<1>(dh3,dh4);
  cm.sew<1>(dh4,dh1);
  write_inc<CM_2>(cm);
  return 0;
}

#include <CGAL/Combinatorial_map.h>

#include <iostream>
#include "inc_graph_utils.h"

typedef CGAL::Combinatorial_map<2> CM_2;
typedef CM_2::Dart_handle Dart_handle;

int main(){
  CM_2 cm;
  Dart_handle dh1 = cm.create_dart();
  Dart_handle dh2 = cm.create_dart();
  Dart_handle dh3 = cm.create_dart();
  cm.sew<1>(dh1,dh2);
  cm.sew<1>(dh2,dh3);
  cm.sew<1>(dh3,dh1);
  export_incidence_graph<CM_2,2>(cm, "test.csv");
  return 0;
}

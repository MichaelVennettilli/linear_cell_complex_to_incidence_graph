#include <CGAL/Combinatorial_map.h>

#include <iostream>
#include "inc_graph_utils.h"

typedef CGAL::Combinatorial_map<3> CM;
typedef CM::Dart_handle Dart_handle;

int main(){
  CM cm;
  cm.make_combinatorial_tetrahedron();
  export_incidence_graph<CM,cm.dimension>(cm, "test.csv");
  return 0;
}

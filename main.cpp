#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

#include "adjacency_incidence_cmap_utils.h"

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
typedef LCC_3::Dart_handle Dart_handle;
typedef LCC_3::Point Point;


int main(){
  LCC_3 lcc;
  lcc.make_tetrahedron(Point(0,0,0),Point(1,0,0),Point(0,1,0),Point(0,0,1));
  export_incidence_graph<LCC_3,lcc.dimension>(lcc, "test_inc.csv");
  export_vertices<LCC_3>(lcc, "test_points.csv");
  return 0;
}

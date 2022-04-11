#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include "minimal_shell_partitioner.h"
#include "adjacency_incidence_cmap_utils.h"

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
typedef LCC_3::Dart_handle Dart_handle;
typedef LCC_3::Point Point;

int main(int argc, char* argv[])
{
  double r_in = 1.0;
  double r_out = 2.0;
  int num_pts = 50;
  if (argc>1){
    num_pts = std::atoi(argv[1]);
  }
  if (argc>2){
		r_in = std::atof(argv[2]);
	}
	if (argc>3){
		r_out = std::atof(argv[3]);
	}

  std::vector<Point> points = random_spherical_points(num_pts);
  LCC_3 shell = generate_shell(points, r_in, r_out);
  lloyd_relaxation(shell, 10, r_in, r_out);
  triangulate_all_faces(shell);
  export_incidence_graph<LCC_3,shell.dimension>(shell, "test_inc.csv");
  export_vertices<LCC_3>(shell, "test_points.csv");
  return EXIT_SUCCESS;
}

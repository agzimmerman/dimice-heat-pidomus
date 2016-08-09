#include <deal.II/grid/tria.h>
namespace Grids {
    using namespace dealii;
    void hemisphere_cylinder_shell(Triangulation<2,2> & tria,
                                 const double inner_radius, // Radius of shell's inside boundary
                                 const double outer_radius, // Radius of shell's outside boundary
                                 const double inner_length, // Length of the inside cylinder
                                 const double outer_length  // Length of the outside cylinder
			         );
}

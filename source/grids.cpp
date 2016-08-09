#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include "grids.h"

namespace Grids {
    using namespace dealii;
    void hemisphere_cylinder_shell(Triangulation<2,2> & tria,
                                 const double inner_radius, // Radius of shell's inside boundary
                                 const double outer_radius, // Radius of shell's outside boundary
                                 const double inner_length, // Length of the inside cylinder
                                 const double outer_length  // Length of the outside cylinder
			         )
    {
        /*
        Based on hyper_cube_with_cylindrical_hole
        Boundary ID's follow the vertex numbering, with 
            0 -> right side of outer spherical boundary
            ... counter-clockwise
            4 -> left side of outer spherical boundary
            5 -> right side of inner spherical boundary
            ... counter-clockwise
            9 -> left side of inner spherical boundary
        Manifold ID's:
            0 -> Hemi-spherical manifold radially centered at the origin
            1 -> Cyldrincal manifoldd
        The origin of the coordinate system is at the radial center of the spherical manifold.
        */
        const int dim = 2;
        Assert(outer_radius > inner_radius, ExcInvalidState());
        Assert(outer_length > inner_length, ExcInvalidState());
        // To ensure valid indexing and to simplify extension to 3D:
        //  - Create an hyper_shell in two dimensions.
        //  - Modify it to match our desired shape.
        GridGenerator::hyper_shell(tria, {0, 0}, inner_radius, outer_radius, 5, true);
        // Rotate the grid ninety degrees to obtain symmetry about the y-axis.
        GridTools::rotate(-numbers::PI/2., tria);
        // Modify half of the grid into a cylindrical shell.
        Triangulation<dim>::active_cell_iterator cell = tria.begin_active(), endc = tria.end();
        std::vector<bool> treated_vertices(tria.n_vertices(), false);
        for (; cell != endc; ++cell) {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
                if (cell->face(f)->at_boundary()) {
                    for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
                        unsigned int vv = cell->face(f)->vertex_index(v);
                        if (treated_vertices[vv] == false) {
                            treated_vertices[vv] = true;
                            switch (vv) {
                                case 1:
                                    cell->face(f)->vertex(v) = Point<dim>(outer_radius,0.);
                                    break;
                                case 2:
                                    cell->face(f)->vertex(v) = Point<dim>(outer_radius,outer_length);
                                    break;
                                case 3:
                                    cell->face(f)->vertex(v) = Point<dim>(-outer_radius,outer_length);
                                    break;
                                case 4:
                                    cell->face(f)->vertex(v) = Point<dim>(-outer_radius,0.);
                                    break;
                                case 6:
                                    cell->face(f)->vertex(v) = Point<dim>(inner_radius,0.);
                                    break;
                                case 7:
                                    cell->face(f)->vertex(v) = Point<dim>(inner_radius,inner_length);
                                    break;
                                case 8:
                                    cell->face(f)->vertex(v) = Point<dim>(-inner_radius,inner_length);
                                    break;
                                case 9:
                                    cell->face(f)->vertex(v) = Point<dim>(-inner_radius,0.);
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // Set boundary ID's corresponding to vertex numbers.
        // This uses the same loop structure as above.
        cell = tria.begin_active();
        std::fill(treated_vertices.begin(), treated_vertices.end(), false);
        double eps = 1e-8*inner_radius;
        for (; cell != endc; ++cell) {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
                if (cell->face(f)->at_boundary()) {
                    double xf = cell->face(f)->center()[0], yf = cell->face(f)->center()[1];
                    for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
                        unsigned int vv = cell->face(f)->vertex_index(v);
                        if (treated_vertices[vv] == false) {
                            switch (vv) {
                                case 0:
                                    if (xf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 1:
                                    if (yf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 2:
                                    if (xf < eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 3:
                                    if (yf < inner_length) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 4:
                                    if (yf < -eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 5:
                                    if (xf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 6:
                                    if (yf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 7:
                                    if (xf < eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 8:
                                    if (yf < inner_length - eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 9:
                                    if (yf < -eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // Label the spherical and cylindrical manifolds.
        cell = tria.begin_active();
        for (; cell!=endc; ++cell) {
            if ((cell->center())[1] < eps) {
                cell->set_all_manifold_ids(0); // Spherical
            }
            else {
                cell->set_all_manifold_ids(1); // Cylindrical
            }
        }
    }
}

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include "my_grids.h"
namespace MyGrids {
    using namespace dealii;
    void make_grid_file() {
        Triangulation<2> grid;
        sphere_cylinder_shell(grid, 0.25, 0.5, 1.0, 1.25);
        std::ofstream out ("grid.vtk");
        GridOut grid_out;
        grid_out.write_vtk (grid, out);
        std::cout << "Grid written to grid.vtk" << std::endl;
        std::ofstream image_out("grid_image.eps");
        grid_out.write_eps(grid, image_out);
    }
    template <int dim, int spacedim>
    void sphere_cylinder_shell(dealii::Triangulation<dim,spacedim> & grid,
                                const double inner_radius, // Radius of shell's inside boundary
                                const double outer_radius, // Radius of shell's outside boundary
                                const double inner_length, // Length of the inside cylinder
                                const double outer_length) // Length of the outside cylinder
                            {
        // Based on create_coarse_grid in http://dealii.org/8.4.1/doxygen/deal.II/step_14.html
        // The origin of the coordinate system is at the center of the spherical manifolds.
        Assert (dim==2, dealii::ExcNotImplemented());
        Assert (outer_radius > inner_radius, dealii::ExcInvalidState());
        Assert (outer_length > inner_length, dealii::ExcInvalidState());
        using dealii::Point;
        static const Point<dim> static_points[] = {
            // Points of the shell's inner boundary:
            Point<dim>(0,-inner_radius),
            Point<dim>(inner_radius,0),
            Point<dim>(inner_radius,inner_length),
            Point<dim>(-inner_radius,inner_length),
            Point<dim>(-inner_radius,0),
            // Points of the shell's outer boundary:
            Point<dim>(0,-outer_radius),
            Point<dim>(outer_radius,0),
            Point<dim>(outer_radius,inner_length),
            Point<dim>(outer_radius,outer_length),
            Point<dim>(inner_radius,outer_length),
            Point<dim>(-inner_radius,outer_length),
            Point<dim>(-outer_radius,outer_length),
            Point<dim>(-outer_radius,inner_length),
            Point<dim>(-outer_radius,0)
        };
        const unsigned int point_count = sizeof(static_points)/sizeof(static_points[0]);
        const std::vector<Point<dim>> points(&static_points[0], &static_points[point_count]);
        const int vertices_per_cell = 4; // ISO C++ forbids variable length array
        const int cell_count = 7;
        // @todo: What are the rules for ordering vertices?
        // It should suffice to have a consistent orientation for each cell, but that does not seem to work.
        // The ordering below arose from tediously guessing and checking.
        // Maybe instead a hyper shell should be created and then modified to our shape. This should guarantee correctness.
        static const int cell_vertices[cell_count][vertices_per_cell] = {
            {4, 13, 0, 5},
            {0, 5, 1, 6},
            {1, 6, 2, 7},
            {2, 7, 9, 8},
            {3, 2, 10, 9},
            {12, 3, 11, 10},
            {13, 4, 12, 3},
        };
        std::vector<dealii::CellData<dim> > cells(cell_count, dealii::CellData<dim>());
        for (unsigned int i=0; i< cell_count; ++i) {
            for (unsigned int j=0; j < vertices_per_cell; ++j) {
                cells[i].vertices[j] = cell_vertices[i][j];
            }
            cells[i].material_id = 0;
        }
        dealii::GridReordering<dim,dim>::reorder_cells(cells, true); // Note: Re-ordering only works with uniformly oriented grid.
        grid.create_triangulation(points,
                                  cells,
                                  dealii::SubCellData());
        // Set boundary id's similar to how it was done in hyper_cube_with_cylindrical_hole.
        double eps = 1e-3*inner_radius;
        auto cell = grid.begin_active(), endc = grid.end();
        for (; cell != endc; ++cell) {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
                if (cell->face(f)->at_boundary()) {
                    dealii::Point<dim> center = cell->face(f)->center();
                    double x = center[0], y = center[1];
                    if ((std::abs(y + inner_radius/2.) < eps) && (std::abs(x) - inner_radius/2.) < eps) {
                        cell->face(f)->set_boundary_id(0); // on the inner spherical boundary
                    }
                    else if ((-eps <= y) & (y <= inner_length + eps) & (std::abs(x) <= inner_radius + eps)) {
                        cell->face(f)->set_boundary_id(1); // on the inner rectangular boundary
                    }
                    else {
                        cell->face(f)->set_boundary_id(2); // on the outer spherical or rectuangular boundary
                    }
                }
            }
        }
        // Set spherical manifolds.
        const dealii::SphericalManifold<dim> manifold_description(dealii::Point<dim>(0,0));
        int spherical_manifold_id = 0;
        grid.set_manifold(spherical_manifold_id, manifold_description);
        cell = grid.begin_active();
        for (; cell!=endc; ++cell) {
            if ((cell->center())[1] < eps) {
                cell->set_all_manifold_ids(spherical_manifold_id);
            }
        }
    }
}

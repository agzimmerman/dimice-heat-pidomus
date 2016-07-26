#ifndef MY_GRID_GENERATOR_H
#define MY_GRID_GENERATOR_H

#include <deal.II/grid/tria.h>

namespace MyGrids
{
    void make_grid_file();
    /**
     * Produce a domain that is the space between two sphere-cylinders.
     * @todo Insert image that shows the geometry and labels the parameters.
     * @note This has only been implemented in 2D.
     */
    template <int dim, int spacedim>
    void sphere_cylinder_shell (dealii::Triangulation<dim,spacedim> &tria,
                                const double inner_radius=0.25,
                                const double outer_radius=0.5,
                                const double inner_length=1.0,
                                const double outer_length=1.25);
}

#endif // MY_GRID_GENERATOR_H

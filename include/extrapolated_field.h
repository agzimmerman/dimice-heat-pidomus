#ifndef _extrapolated_field_h_
#define _extrapolated_field_h_

#include <deal.II/numerics/fe_field_function.h>
/**
 * @brief Extends the FEFieldFunction class with nearest neighbor extrapolation.
 *
 * @detail
 *
 *    This is perhaps useful for interpolating the function onto a similar domain
 *    with non-conforming boundaries, e.g. if the domain has only been slightly shifted or rotated.
 *
 * @author Alexander Zimmerman 2016
*/
namespace MyFunctions
{
template< int dim,
          typename DoFHandlerType=DoFHandler<dim>,
          typename VectorType=Vector<double> >
class ExtrapolatedField : public Functions::FEFieldFunction<dim,DoFHandlerType,VectorType> {
    public:
	ExtrapolatedField(const DoFHandlerType &_dof_handler, const VectorType &f)
	    : Functions::FEFieldFunction<dim,DoFHandlerType,VectorType> (_dof_handler, f),
            dof_handler(&_dof_handler, "ExtrapolatedField") {}
	virtual double extrapolated_value(const Point<dim> &p,
                                          const unsigned int component = 0);
    private:
        SmartPointer<const DoFHandlerType,ExtrapolatedField<dim,DoFHandlerType,VectorType>>	  
	    dof_handler;
	Point<dim> get_nearest_boundary_vertex(const Point<dim> &point);
};  

template<int dim, typename DoFHandlerType, typename VectorType>
double ExtrapolatedField<dim,DoFHandlerType,VectorType>::
	extrapolated_value(const Point<dim> &point,
                           const unsigned int) {
    double val;
    try {
        val = this->value(point);
    }
    catch (VectorTools::ExcPointNotAvailableHere) {
	val = this->value(this->get_nearest_boundary_vertex(point));
    }
    return val;
}

template <int dim, typename DoFHandlerType, typename VectorType>
Point<dim> ExtrapolatedField<dim,DoFHandlerType,VectorType>::
	get_nearest_boundary_vertex(const Point<dim> &point) {
    double arbitrarily_large_number = 1.e32;
    double nearest_distance = arbitrarily_large_number;
    Point<dim> nearest_vertex;
    for (auto cell : dof_handler->active_cell_iterators()) {
	if (!cell->at_boundary()) {
	    continue;
	}
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
	    if (!cell->face(f)->at_boundary()) {
		continue;
	    }
            for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
		Point<dim> vertex = cell->face(f)->vertex(v);
		double distance = (point - vertex).norm_square();
		if (distance < nearest_distance) {
		    nearest_vertex = vertex;
		    nearest_distance = distance;		
		}
	    }
	}
    }
    return nearest_vertex;
}
}

#endif


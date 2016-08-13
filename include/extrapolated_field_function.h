#include <deal.II/base/function>
/*
 * @brief Extends the FEFieldFunction class with nearest neighbor extrapolation outside of the domain.
 *
 * @detail
 *
 *    This is perhaps useful for interpolating the function onto a similar domain
 *    with non-conforming boundaries, e.g. if the domain has only been slightly shifted or rotated.
 *
 * @author Alexander Zimmerman 2016
*/
namespace Functions
{
template<int dim,
         typename DoFHandlerType=DoFHanlder<dim>,
         typename VectorType=Vector<double>>
class ExtrapolatedField : public Function<dim> {
    public:
	ExtrapolatedField() : Function<dim>() {};
	ExtrapolatedField(&DoF,
		          &Vector);
	virtual double value(const Point<dim> &p, 
			     const unsigned int component = 0) const;
	void set_field(&DoF,
		       &Vector);
    private:
	Functions::FEFieldFunction<dim,DoF,Vector> field;
	Point<dim> get_nearest_boundary_vertex(const DoF &dof_handler,
                                               const Point<dim> &point);
};

template<int dim, typename Vector>
ExtrapolatedField::ExtrapolatedField(DoFHandler<dim,dim> &dof_handler,
                                     Vector f) {
    ExtrapolatedField();
    set_field(dof_handler, f);
}

template<int dim, DoFHandler<dim,dim>, typename Vector>
void ExtrapolatedField<dim>::set_field(DoFHandler<dim,dim> &dof_handler,
                                       Vector &f) {
    auto field(dof_handler, f);
}

template<int dim>
double ExtrapolatedField<dim>::value (const Point<dim> &point,
                                      const unsigned int) const {
    double val;
    try {
        val = field.value(point)
    }
    catch {
	val = field.value(get_nearest_boundary_vertex(field.dof_handler,
                                                      point));
    }
    return val;
}

template <int dim>
Point<dim> ExtrapolatedField::get_nearest_boundary_vertex(const DoFHandler<dim,dim> &dof_handler,
                                                          const Point<dim> &point) {
    double arbitrarily_large_number = 1.e32;
    double nearest_distance = arbitrarily_large_number;
    Point<dim> nearest_vertex;
    for (auto cell : dof_handler.active_cell_iterators()) {
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
	    // Check only the vertices that are on boundary face
	    if (!cell->face(f)->at_boundary()) {
		continue;
	    }
            for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
		Point<dim> vertex = cell->face(f)->vertex(v)
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

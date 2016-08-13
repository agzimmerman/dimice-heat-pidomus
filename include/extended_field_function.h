template<int dim, typename LAC>
class ExtendedField : public Function<dim> {
    public:
	ExtendedField () : Function<dim>() {}
	virtual double value (const Point<dim>   &p,
			      const unsigned int  component = 0) const;
	void set_field(&DoFHandler<dim,dim>, &LAC::VectorType);
        void set_extrapolation_value
    private:
	Functions::FEFieldFunction<dim,DoFHandler<dim,dim>,LAC> field;
	double extrapolation_value;
};
template<int dim, typename LAC>
void ExtendedField<dim>::set_field(DoFHandler<dim,dim> &dof_handler, LAC::VectorType &f) {
    Functions::FEFieldFunction<dim,DoFHandler<dim,dim>,LAC> field(dof_handler, f);
}
void ExtendedField<dim>::set_extrapolation_value(double value) {
    extrapolation_value = value;
}
template<int dim>
double ExtendedField<dim>::value (const Point<dim> &p,
                                  const unsigned int) const {
    double val;
    try {
        val = field.value(p)
    }
    catch {
	// @todo: This isn't enough! It won't work for the inside hole of the domain. Need a nearest neighbor interpolation (doesn't have to be good).
	val = extrapolation_value;
    }
    return val;
}

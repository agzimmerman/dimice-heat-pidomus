#include <pidomus.h>
#include <interfaces/poisson_problem.h>
#include <string>
#include <heat_problem.hpp>
template <int dim> HeatProblem<dim>::HeatProblem() {
    _dim = dim;
    PoissonProblem<dim,dim,LADealII> _pde_interface;
}
template <int dim> void HeatProblem<dim>::solve(std::string parameter_file_path) {
    piDoMUS<dim,dim,LADealII> solver("pi-DoMUS", this->_pde_interface);
    ParameterAcceptor::initialize(parameter_file_path, "used_parameters.prm");
    solver.run();
}

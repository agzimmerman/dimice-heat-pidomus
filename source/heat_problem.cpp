#include <pidomus.h>
#include <interfaces/poisson_problem.h>
#include <string>
#include <heat_problem.hpp>
HeatProblem::HeatProblem() {
    PoissonProblem<this->_dim,this->_dim,LADealII> _pde_interface;
}
void HeatProblem::solve(std::string parameter_file_path) {
    piDoMUS<this->_dim,this->_dim,LADealII> solver("pi-DoMUS", this->_pde_interface);
    ParameterAcceptor::initialize(parameter_file_path, "used_parameters.prm");
    solver.run();
}

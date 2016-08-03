#include <heat_problem.hpp>
HeatProblem::HeatProblem(std::string parameter_file_path) {
    _parameter_file_path = parameter_file_path;
    PoissonProblem<this->_dim,this->_dim,LADealII> _pde_interface;
}
void HeatProblem::solve() {
    piDoMUS<this->_dim,this->_dim,LADealII> solver("pi-DoMUS", this->_pde_interface);
    ParameterAcceptor::initialize(_parameter_file_path, "used_parameters.prm");
    solver.run();
}

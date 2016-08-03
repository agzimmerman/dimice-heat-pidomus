#include <heat_problem.hpp>
#define _unused(x) ((void)(x)) // http://stackoverflow.com/questions/777261/avoiding-unused-variables-warnings-when-using-assert-in-a-release-build
HeatProblem::HeatProblem(const std::string parameter_file_path, const int mpi_size) {
    // Initialize MPI
    // This approach emulates argc and argv to use the deal.II utility MPI_InitFinalize
    int argc = 3;
    std::vector<char*> cstrings;
    std::string my_name = "HeatProblem";
    cstrings.push_back(const_cast<char*>(my_name.c_str()));
    cstrings.push_back(const_cast<char*>("-n"));
    cstrings.push_back(const_cast<char*>((std::to_string(mpi_size)).c_str()));
    char** argv = cstrings.data();
    Utilities::MPI::MPI_InitFinalize _mpi_init(argc, argv, numbers::invalid_unsigned_int);
    // Initialize the PDE problem.
    _parameter_file_path = parameter_file_path;
    PoissonProblem<this->_dim,this->_dim,LADealII> _pde_interface;
}
void HeatProblem::solve() {
    piDoMUS<this->_dim,this->_dim,LADealII> solver("pi-DoMUS", this->_pde_interface);
    ParameterAcceptor::initialize(_parameter_file_path, "used_parameters.prm");
    solver.run();
}

#include <pidomus.h>
#include "interfaces/poisson_problem_signals.h"
#define _unused(x) ((void)(x)) // http://stackoverflow.com/questions/777261/avoiding-unused-variables-warnings-when-using-assert-in-a-release-build
int main (int argc, char *argv[]) {
    const int dim = 2;
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, numbers::invalid_unsigned_int);
    std::cout << argv[0] << std::endl;
    std::string parameter_file_path = "";
    if (argc > 1) {
        parameter_file_path = argv[1];
        std::cout << "Using parameter input file at " << parameter_file_path << std::endl;
    }
    else {
        std::string exe_path = argv[0];
        parameter_file_path = exe_path.substr(0, exe_path.find_last_of("\\/"))+"/../inputs/sphere.prm";
        std::cout << "Using default parameters at " << parameter_file_path << std::endl;
    }
    PoissonProblem<dim,dim,LADealII> p;
    piDoMUS<dim,dim,LADealII> solver("pi-DoMUS", p);
    ParameterAcceptor::initialize(parameter_file_path, "used_parameters.prm");
    solver.run();
    _unused(mpi_init);
    return 0;
}

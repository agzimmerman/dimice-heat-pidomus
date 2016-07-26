#include <pidomus.h>
#include <interfaces/poisson_problem.h>
#include "my_grids.h"
#define _unused(x) ((void)(x)) // http://stackoverflow.com/questions/777261/avoiding-unused-variables-warnings-when-using-assert-in-a-release-build
int main (int argc, char *argv[]) {
    const int dim = 2;
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, numbers::invalid_unsigned_int);
    MyGrids::make_grid_file();
    std::string parameter_file_path = "";
    if (argc > 1) {
        parameter_file_path = argv[1];
        std::cout << "Using parameter input file at " << parameter_file_path << std::endl;
    }
    else {
        std::cout << "Using default parameters." << std::endl;
    }
    PoissonProblem<dim,dim,LADealII> p;
    piDoMUS<dim,dim,LADealII> solver("pi-DoMUS", p);
    ParameterAcceptor::initialize(parameter_file_path, "used_parameters.prm");
    solver.run();
    _unused(mpi_init);
    return 0;
}
#include <iostream>
#include <heat_problem.hpp>
#define _unused(x) ((void)(x)) // http://stackoverflow.com/questions/777261/avoiding-unused-variables-warnings-when-using-assert-in-a-release-build
int main (int argc, char *argv[]) {
    // Parse command line inputs for parameter file path.
    std::string parameter_file_path = "";
    if (argc > 1) {
        parameter_file_path = argv[1];
        std::cout << "Using parameter input file at " << parameter_file_path << std::endl;
    }
    else {
        std::string exe_path = argv[0];
        parameter_file_path = exe_path.substr(0, exe_path.find_last_of("\\/"))+"/../inputs/sphere-cylinder.prm";
        std::cout << "Using default parameters at " << parameter_file_path << std::endl;
    }
    // Run the heat problem for the given parameter file.
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, numbers::invalid_unsigned_int);
    HeatProblem<2> heat;
    heat.solve(parameter_file_path);
    _unused(mpi_init);
    return 0;
}

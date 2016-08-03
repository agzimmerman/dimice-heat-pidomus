#include <heat_problem.hpp>
#define _unused(x) ((void)(x)) // http://stackoverflow.com/questions/777261/avoiding-unused-variables-warnings-when-using-assert-in-a-release-build
HeatProblem::HeatProblem(std::vector<std::string> args) {
    // Emulate argc and argv so that we can use them to initiliaze MPI.
    std::vector<char*> cstrings;
    int argc = args.size();
    for(size_t i = 0; i < argc; ++i) {
        cstrings.push_back(const_cast<char*>(args[i].c_str()));
    }
    char** argv = cstrings.data();
    // Initialize MPI
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, numbers::invalid_unsigned_int);
    // Parse command line inputs for parameter file path.
    std::string parameter_file_path = "";
    if (argc > 1) {
        parameter_file_path = argv[1];
        std::cout << "Using parameter input file at " << parameter_file_path << std::endl;
    }
    else {
        std::string exe_path = argv[0];
        parameter_file_path = exe_path.substr(0, exe_path.find_last_of("\\/"))+"/../inputs/sphere-cylinder.prm";
    }
    _parameter_file_path = parameter_file_path;
    PoissonProblem<this->_dim,this->_dim,LADealII> _pde_interface;
    _unused(mpi_init);
}
void HeatProblem::solve() {
    piDoMUS<this->_dim,this->_dim,LADealII> solver("pi-DoMUS", this->_pde_interface);
    ParameterAcceptor::initialize(_parameter_file_path, "used_parameters.prm");
    solver.run();
}

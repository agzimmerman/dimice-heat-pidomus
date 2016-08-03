#include <iostream>
#include <heat_problem.hpp>

int main (int argc, char *argv[]) {
std::cout << "Using default parameters at " << parameter_file_path << std::endl;
    }
    // Run the heat problem for the given parameter file.
    HeatProblem heat(parameter_file_path);
    heat.solve();
    return 0;
}

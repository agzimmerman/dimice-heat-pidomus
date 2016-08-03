#ifndef HEAT_PROBLEM_H
#define HEAT_PROBLEM_H
#include <string>
#include <pidomus.h>
#include <interfaces/poisson_problem.h>
class HeatProblem {
public:
    HeatProblem(
        const std::string parameter_file_path="/../inputs/sphere-cylinder.prm",
        const int mpi_size=1);
    void solve();
private:
    Utilities::MPI::MPI_InitFinalize _mpi_init;
    static const int _dim = 2;
    std::string _parameter_file_path;
    PoissonProblem<_dim,_dim,LADealII>_pde_interface;
};
#endif // HEAT_PROBLEM_H

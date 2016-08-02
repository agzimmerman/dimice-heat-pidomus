#ifndef HEAT_PROBLEM_H
#define HEAT_PROBLEM_H
#include <interfaces/poisson_problem.h>
#include <string>
template <int dim> class HeatProblem {
public:
    HeatProblem();
    void solve(std::string parameter_file_path="");
private:
    const int _dim;
    PoissonProblem<dim,dim,LADealII>_pde_interface;
};
#endif // HEAT_PROBLEM_H

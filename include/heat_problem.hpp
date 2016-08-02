#ifndef HEAT_PROBLEM_H
#define HEAT_PROBLEM_H
#include <interfaces/poisson_problem.h>
#include <string>
class HeatProblem {
public:
    HeatProblem();
    void solve(std::string parameter_file_path="");
private:
    static const int _dim = 2;
    PoissonProblem<_dim,_dim,LADealII>_pde_interface;
};
#endif // HEAT_PROBLEM_H

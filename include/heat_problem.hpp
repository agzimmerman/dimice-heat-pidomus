#ifndef HEAT_PROBLEM_H
#define HEAT_PROBLEM_H
#include <interfaces/poisson_problem.h>
#include <string>
class HeatProblem {
public:
    HeatProblem(std::string parameter_file_path="");
    void solve();
private:
    static const int _dim = 2;
    std::string _parameter_file_path;
    PoissonProblem<_dim,_dim,LADealII>_pde_interface;
};
#endif // HEAT_PROBLEM_H

#include "pidomus.h"
#include "interfaces/poisson_problem.h"
#include <string>
using namespace dealii;
int main (int argc, char *argv[])
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
                                                      numbers::invalid_unsigned_int);
  const int dim = 2;
  const int spacedim = 3;

  std::string parameter_filepath = "";
  if (argc > 1) {
    parameter_filepath = argv[1];
	std::cout << "Using parameter input file at " << parameter_filepath << std::endl;
  }
  else {
    std::cout << "Using default parameters." << std::endl;
  }
  PoissonProblem<dim,spacedim,LADealII> p;
  piDoMUS<dim,spacedim,LADealII> solver ("pidomus",p);
  ParameterAcceptor::initialize(parameter_filepath, "used_parameters.prm");


  solver.run ();

  auto sol = solver.get_solution();
  for (unsigned int i = 0; i<sol.size(); ++i)
    deallog << sol[i] << std::endl;

  return 0;
}

#ifndef _pidoums_poisson_h_
#define _pidoums_poisson_h_

#include <iostream>
#include <fstream>
#include <string>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include "pde_system_interface.h"


using namespace dealii;

template <int dim, int spacedim, typename LAC=LATrilinos>
class PoissonProblem : public PDESystemInterface<dim,spacedim, PoissonProblem<dim,spacedim,LAC>, LAC>
{

public:
  ~PoissonProblem () {};
  PoissonProblem ();

  // interface with the PDESystemInterface :)


  template <typename EnergyType, typename ResidualType>
  void energies_and_residuals(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
                              FEValuesCache<dim,spacedim> &scratch,
                              std::vector<EnergyType> &energies,
                              std::vector<std::vector<ResidualType> > &local_residuals,
                              bool compute_only_system_terms) const;


  void compute_system_operators(const std::vector<shared_ptr<LATrilinos::BlockMatrix> >,
                                LinearOperator<LATrilinos::VectorType> &,
                                LinearOperator<LATrilinos::VectorType> &,
                                LinearOperator<LATrilinos::VectorType> &) const;


  virtual void connect_to_signals() const
  {
    // first of all we get the struct Signals from pidomus
    auto &signals = this->get_signals();

    signals.postprocess_newly_created_triangulation.connect([&]
	    (typename parallel::distributed::Triangulation<dim,spacedim> &tria) {
	tria.save(serial_coarse_tria_file_name.c_str());
    });
        

    // Connect to signal to serialize data before returning from pi-DoMUS.
    signals.use_solution_before_return.connect([&](parallel::distributed::Triangulation<dim,spacedim> &tria,
						   DoFHandler<dim,spacedim> &dof_handler,
                                                   typename LAC::VectorType &solution,
                                                   typename LAC::VectorType &solution_dot) {
	std::cout << "Connected to signals.serialize_before_return" << std::endl;
	parallel::distributed::SolutionTransfer<dim,typename LAC::VectorType> sol_trans(dof_handler);
	sol_trans.prepare_serialization(solution);
	tria.save("serialized_solution.txt");
	sol_trans.prepare_serialization(solution_dot);
	tria.save("serialized_solution_dot.txt");
    });

    // Connect to signal to modify initial conditions.
    signals.deserialize_initial_conditions.connect([&](DoFHandler<dim,dim> &dof_handler,
					       typename LAC::VectorType &solution,
                                               typename LAC::VectorType &solution_dot) {
	std::cout << "Connected to signals.deserialize_initial_conditions" << std::endl;
	// If a serialized solution exists, then initialize with that solution.
	{	
	    std::ifstream file_to_check("serialized_solution.txt");
	    if (!file_to_check.good()) {
	        return;
	    }
        }
	std::cout << "Initializing from serialized solution." << std::endl;
        // Load serialized data
	double r0 = 0.25, r1 = 0.5, l0 = 1.0, l1 = 1.25;
        Point<3> trans;
	trans[0] = 0;
	trans[1] = 0;
	trans[2] = 0;
	parallel::distributed::Triangulation<dim,spacedim> parallel_coarse_tria(MPI_COMM_WORLD);
        GridGenerator::hemisphere_cylinder_shell(parallel_coarse_tria, r0, r1, l0, l1, trans);
        parallel_coarse_tria.load("serialized_solution.txt");
	parallel::distributed::SolutionTransfer<dim,typename LAC::VectorType> sol_trans(dof_handler);
	sol_trans.deserialize(solution);
        parallel_coarse_tria.load("serialized_solution_dot.txt");
	sol_trans.deserialize(solution_dot);
    });

  }

private:
  mutable shared_ptr<TrilinosWrappers::PreconditionJacobi> preconditioner;
  std::string serial_coarse_tria_file_name = "serialized_coarse_triangulation.txt";
};

template <int dim, int spacedim, typename LAC>
PoissonProblem<dim,spacedim, LAC>::
PoissonProblem():
  PDESystemInterface<dim,spacedim,PoissonProblem<dim,spacedim,LAC>, LAC >("Poisson problem",
      1,1,
      "FESystem[FE_Q(1)]",
      "u","1")
{}



template <int dim, int spacedim, typename LAC>
template <typename EnergyType, typename ResidualType>
void
PoissonProblem<dim,spacedim,LAC>::
energies_and_residuals(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
                       FEValuesCache<dim,spacedim> &fe_cache,
                       std::vector<EnergyType> &,
                       std::vector<std::vector<ResidualType> > &local_residuals,
                       bool compute_only_system_terms) const
{

  const FEValuesExtractors::Scalar s(0);

  ResidualType rt = 0; // dummy number to define the type of variables
  this->reinit (rt, cell, fe_cache);
  auto &uts = fe_cache.get_values("solution_dot", "u", s, rt);
  auto &gradus = fe_cache.get_gradients("solution", "u", s, rt);

  const unsigned int n_q_points = uts.size();
  auto &JxW = fe_cache.get_JxW_values();

  auto &fev = fe_cache.get_current_fe_values();

  for (unsigned int q=0; q<n_q_points; ++q)
    {
      auto &ut = uts[q];
      auto &gradu = gradus[q];
      for (unsigned int i=0; i<local_residuals[0].size(); ++i)
        {
          auto v = fev[s].value(i,q);
          auto gradv = fev[s].gradient(i,q);
          local_residuals[0][i] += (
                                     ut*v
                                     +
                                     gradu*gradv
                                   )*JxW[q];
        }

      (void)compute_only_system_terms;

    }

}


template <int dim, int spacedim, typename LAC>
void
PoissonProblem<dim,spacedim,LAC>::compute_system_operators(const std::vector<shared_ptr<LATrilinos::BlockMatrix> > matrices,
                                                           LinearOperator<LATrilinos::VectorType> &system_op,
                                                           LinearOperator<LATrilinos::VectorType> &prec_op,
                                                           LinearOperator<LATrilinos::VectorType> &) const
{

  preconditioner.reset  (new TrilinosWrappers::PreconditionJacobi());
  preconditioner->initialize(matrices[0]->block(0,0));

  auto A  = linear_operator<LATrilinos::VectorType::BlockType>( matrices[0]->block(0,0) );

  LinearOperator<LATrilinos::VectorType::BlockType> P_inv;

  P_inv = linear_operator<LATrilinos::VectorType::BlockType>(matrices[0]->block(0,0), *preconditioner);

  auto P00 = P_inv;

  // ASSEMBLE THE PROBLEM:
  system_op  = block_operator<1, 1, LATrilinos::VectorType>({{
      {{ A }}
    }
  });

  prec_op = block_operator<1, 1, LATrilinos::VectorType>({{
      {{ P00}} ,
    }
  });
}

#endif

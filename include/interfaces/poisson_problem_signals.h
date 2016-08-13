#ifndef _pidoums_poisson_h_
#define _pidoums_poisson_h_

#include <iostream>
#include <fstream>
#include <string>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <deal.II/numerics/fe_field_function.h>

#include "pde_system_interface.h"
#include "grids.h"
#include "extrapolated_field_function.h"


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

    // Connect to signal to use custom grid function.
    signals.postprocess_newly_created_triangulation.connect(
            [&](typename parallel::distributed::Triangulation<dim,spacedim> &tria) {
	double inner_radius = 0.25, outer_radius = 0.5, inner_length = 1.0, outer_length = 1.25;
	tria.clear();
	Grids::hemisphere_cylinder_shell(tria,
	                                 inner_radius,
                                         outer_radius,
                                         inner_length,
                                         outer_length);
        // Attach spherical manifold
	static const SphericalManifold<2> manifold_description(Point<2>(0,0));
	tria.set_manifold(0, manifold_description);
    });
    // Connect to signal to serialize data before returning from pi-DoMUS.
    signals.serialize_before_return.connect([&](DoFHandler<dim,dim> &dof_handler,
                                                LATrilinos::VectorType &solution,
                                                LATrilinos::VectorType &solution_dot) {
	{
	    std::ofstream fs("serialized_dof_handler.txt");
            boost::archive::text_oarchive archive(fs);
            archive << dof_handler;
	}
        {
	    std::ofstream fs("serialized_solution.txt");
            boost::archive::text_oarchive archive(fs);
            archive << solution;
	}
	{
	    std::ofstream fs("serialized_solution_dot.txt");
            boost::archive::text_oarchive archive(fs);
            archive << solution_dot;
	}
    });
    // Connect to signal to modify initial conditions.
    signals.fix_initial_conditions.connect([&](DoFHandler<dim,dim> &dof_handler,
					       LATrilinos::VectorType &solution,
                                               LATrilinos::VectorType &solution_dot) {
	// If a serialized solution exists, then initialize with that solution.
	{	
	    std::ifstream file_to_check("serialized_solution.txt");
	    if (!file_to_check.good()) {
	        return;
	    }
        }
        // Load serialized data
	DoFHandler<dim,dim> old_dof_handler;
	LAC old_solution, old_solution_dot;
	{
	    std::ifstream fs("serialized_dof_handler.txt");
            boost::archive::text_iarchive archive(fs);
            archive >> old_dof_handler;
	}
	{
	    std::ifstream fs("serialized_solution.txt");
            boost::archive::text_iarchive archive(fs);
            archive >> old_solution;
	}
	{
	    std::ifstream fs("serialized_solution_dot.txt");
            boost::archive::text_iarchive archive(fs);
            archive >> old_solution_dot;
	}
	// Transform the serialized domain per the new state vector.

	// Make the FE field functions	
        ExtrapolatedField<dim,LAC> solution_field(old_dof_handler, old_solution);
	ExtrapolatedField<dim,LAC> solution_dot_field(old_dof_handler, old_solution_dot);
	// VectorTools::interpolate transformed old solution onto current FE space
	VectorTools::interpolate(dof_handler, solution_field, solution);
   	VectorTools::interpolate(dof_handler, solution_dot_field, solution_dot);
    });
  }
private:
  mutable shared_ptr<TrilinosWrappers::PreconditionJacobi> preconditioner;
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

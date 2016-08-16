#ifndef _pidoums_poisson_h_
#define _pidoums_poisson_h_

#include <iostream>
#include <fstream>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include "extrapolated_field.h"
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
	Triangulation<dim,spacedim> serial_tria;
	serial_tria.copy_triangulation(tria);
	std::ofstream file_stream(serial_coarse_tria_file_name, std::ios::out);
	if (!file_stream.good()) {
	    throw std::runtime_error("Error while opening the file: " + serial_coarse_tria_file_name);
	}
	boost::archive::text_oarchive output_archive(file_stream);
	output_archive << serial_tria;
    });
        

    // Connect to signal to serialize data before returning from pi-DoMUS.
    signals.use_solution_before_return.connect([&](const MPI_Comm &comm,
				  	           parallel::distributed::Triangulation<dim,spacedim> &tria,
						   DoFHandler<dim,spacedim> &dof_handler,
                                                   typename LAC::VectorType &solution,
                                                   typename LAC::VectorType &solution_dot) {
	std::cout << "Connected to signals.serialize_before_return" << std::endl;
	// Following example in save method at
	// https://github.com/ORNL-CEES/Cap/blob/master/cpp/source/deal.II/supercapacitor.templates.h#L415
	unsigned int const n_blocks = solution.n_blocks();
	dealii::IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
	dealii::IndexSet locally_relevant_dofs;
	dealii::DoFTools::extract_locally_relevant_dofs(dof_handler,
							locally_relevant_dofs);
	std::vector<dealii::IndexSet> locally_owned_index_sets(n_blocks,
							       locally_owned_dofs);
	std::vector<dealii::IndexSet> locally_relevant_index_sets(n_blocks,
								  locally_relevant_dofs);
	typename LAC::VectorType ghosted_solution(locally_owned_index_sets,
					 locally_relevant_index_sets,
					 comm);
	ghosted_solution = solution;
	dealii::parallel::distributed::SolutionTransfer<dim,typename LAC::VectorType>
		solution_transfer(dof_handler);
	solution_transfer.prepare_serialization(ghosted_solution);
//	solution_transfer.prepare_serialization(solution);
	tria.save(solution_file_name.c_str());
	ghosted_solution = solution_dot;
	solution_transfer.prepare_serialization(ghosted_solution);
//	solution_transfer.prepare_serialization(solution);
	tria.save(solution_dot_file_name.c_str());
	// I think I also need to save this dof_handler for my problem.
	std::ofstream file_stream(dof_handler_file_name, std::ios::out);
	if (!file_stream.good()) {
	    throw std::runtime_error("Error while opening the file: " + dof_handler_file_name);
	}
	boost::archive::text_oarchive output_archive(file_stream);
	output_archive << dof_handler;
    });

    // Connect to signal to modify initial conditions.
    signals.deserialize_initial_conditions.connect([&](const MPI_Comm &comm,
						       FiniteElement<dim,dim> &fe,
						       DoFHandler<dim,dim> &dof_handler,
					               typename LAC::VectorType &solution,
                                                       typename LAC::VectorType &solution_dot) {
	std::cout << "Connected to signals.deserialize_initial_conditions" << std::endl;
	// If a serialized solution exists, then initialize with that solution.
	{	
	    std::ifstream file_to_check("serialized_solution");
	    if (!file_to_check.good()) {
	        return;
	    }
        }
	std::cout << "Initializing from serialized solution." << std::endl;
	// Following example in save method at
	// https://github.com/ORNL-CEES/Cap/blob/master/cpp/source/deal.II/supercapacitor.templates.h#L415
        // Load serialized data
	// Following instructions from Bruno Turcksin's post on the mailing list:
	// 1) Serialize the parallel::Triangulation of the coarse mesh (already did this)
	// 2) Do some refinements, and then use save on the parallel::Triangulation (already did this)
	// To load the data
	// 1) Deserialize the coarse mesh using a Triangulation (NOT a parallel:Triangulation)
	//if (boost::filesystem::exists(serial_coarse_tria_file_name) == false) {
	//    throw std::runtime_error("The file " + serial_coarse_tria_file_name + " does not exist.");
	//}
//	std::ifstream input_file_stream(serial_coarse_tria_file_name, std::ios::binary);
	std::ifstream input_file_stream(serial_coarse_tria_file_name, std::ios::in);
	if (!input_file_stream.good()) {
	    throw std::runtime_error("Error while opening the file: " + serial_coarse_tria_file_name);
	}
//	boost::iostreams::filtering_streambuf<boost::iostreams::input> compressed_in;
//	compressed_in.push(boost::iostreams::zlib_decompressor());
//	compressed_in.push(input_file_stream);
//	boost::archive::binary_iarchive input_archive(compressed_in);
	boost::archive::text_iarchive input_archive(input_file_stream);
        Triangulation<dim> serial_old_coarse_tria;
	input_archive >> serial_old_coarse_tria;
	// 2) Use copy_triangulation to copy the Triangulation to the parallel_Triangulation
        parallel::distributed::Triangulation<dim,spacedim> old_tria(comm);
	old_tria.copy_triangulation(serial_old_coarse_tria);
	// 3) Use load on parallel::Triangulation to get the mesh correctly refined
	//if (boost::filesystem::exists(solution_file_name) == false) {
	//    throw std::runtime_error("The file " + solution_file_name + " does not exist.");
	//}
	old_tria.load(solution_file_name.c_str());
	// That's the end of Bruno's instructions.
	//
	// Now we should be able to use the SolutionTransfer class
//	std::ifstream dof_input_file_stream(dof_handler_file_name, std::ios::in);
//	if (!dof_input_file_stream.good()) {
//	    throw std::runtime_error("Error while opening the file: " + dof_handler_file_name);
//	}
//	boost::archive::text_iarchive dof_input_archive(dof_input_file_stream);
	// Maybe I have to distribute_dofs first: https://groups.google.com/forum/#!searchin/dealii/save$20DoFHandler|sort:relevance/dealii/OI73AWfrN5w/PEBrx_3FICkJ
	// http://dealii.org/8.4.1/doxygen/deal.II/classDoFHandler.html#a553ca864aaf70330d9be86bc78f36d1e
//	shared_ptr<DoFHandler<dim,spacedim>> old_dof_handler;
//	old_dof_handler = SP(new DoFHandler<dim,spacedim>(old_tria));
//	old_dof_handler->distribute_dofs(fe);
//	dof_input_archive >> *old_dof_handler;
	// Or maybe set up dof this way instead: https://groups.google.com/forum/#!searchin/dealii/load$20DoFHandler|sort:relevance/dealii/eWyFll9gp4k/6X-90iwRJAAJ
	dealii::DoFHandler<dim,spacedim> old_dof_handler(old_tria);
	old_dof_handler.distribute_dofs(fe);
	std::ifstream dof_input_file_stream(dof_handler_file_name, std::ios::in);
	if (!dof_input_file_stream.good()) {
	    throw std::runtime_error("Error while opening the file: " + dof_handler_file_name);
	}
	boost::archive::text_iarchive dof_input_archive(dof_input_file_stream);
	dof_input_archive >> old_dof_handler;
	dealii::parallel::distributed::SolutionTransfer<dim,typename LAC::VectorType>
		sol_trans(old_dof_handler);
	//
//	parallel::distributed::SolutionTransfer<dim,typename LAC::VectorType> sol_trans(*old_dof_handler);
	unsigned int ndof = old_dof_handler.n_dofs();
	typename LAC::VectorType old_solution(ndof), old_solution_dot(ndof);
	sol_trans.deserialize(old_solution);
	// Repeat for solution_dot
	//if (boost::filesystem::exists(solution_dot_file_name) == false) {
	//    throw std::runtime_error("The file " + solution_dot_file_name + " does not exist.");
	//}
	old_tria.load(solution_dot_file_name.c_str());
	sol_trans.deserialize(old_solution_dot);
	// Make the FE field functions	
//      MyFunctions::ExtrapolatedField<dim,DoFHandler<dim,dim>,typename LAC::VectorType>
//	   solution_field(*old_dof_handler, old_solution);
//	MyFunctions::ExtrapolatedField<dim,DoFHandler<dim,dim>,typename LAC::VectorType>
//	    solution_dot_field(*old_dof_handler, old_solution_dot);
	// Interpolate transformed old solution onto current FE space
//	VectorTools::interpolate(dof_handler, solution_field, solution);
//   	VectorTools::interpolate(dof_handler, solution_dot_field, solution_dot);
    });

  }

private:
  mutable shared_ptr<TrilinosWrappers::PreconditionJacobi> preconditioner;
  std::string const serial_coarse_tria_file_name = "serialized_coarse_triangulation";
  std::string const solution_file_name = "serialized_solution";
  std::string const solution_dot_file_name = "serialized_solution_dot";
  std::string const dof_handler_file_name = "dof_handler";
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

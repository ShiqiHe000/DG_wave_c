// Doxygen mainpage
/**
 * \author Shiqi He
 * \version 1.0
 * \mainpage Quardtree AMR Discontinuous Galerkin Wave equation
 * \section intro_sec Introduction 
 * A hp-adaptive parallel Discontinuous Galerkin Spectral Element Solver
 * for wave equation. The solver uses quardtree data structure to 
 * achieve AMR. Space-filling curve (HIlbert curve) is invoked for 
 * dynamic load-balancing. 
 * \section compile_sec Compilation
 * Discribe how to compile the code in the future. 
 * \subsection prerequiste Prerequisite
 * e.g. g++ compiler
 * \subsection cmake cmake
 * use cmake
 */
 
#include <iostream>
#include <mpi.h>
#include "dg_mpi.h"
#include "dg_param.h"
#include "dg_prepare_hilbert_scheme.h"
#include "dg_start_parallel.h"
#include "dg_advection_diffusion_driver.h"
#include "dg_basis.h"	// for testing
#include "dg_constructor.h" 	// testing
#include "dg_poly_level_and_order.h"	// test
#include "dg_nodal_2d_storage.h"	//test
#include "dg_single_index.h"	//test
int main(int argc, char *argv[]){
	
	char** a = argv;	
	
	// initialize mpi
	Start_mpi(&argc, a);
	
	// prepare Hilbert curve
	Hilber_numbering();

	// start parallel process
	Start_parallel();

	// start the game
	Driver_for_DG_approximation();
	// testing
//	Construct_basis();
//	int level_max = Poly_order_to_level(grid::nmin, grid::nmax);
//	for(int i = 0; i < level_max + 1; ++i){
////		
//		int porder = Poly_level_to_order(grid::nmin, i);
////		
//		for(int j = 0; j <= porder; ++j ){
////			for(int k = 0; k <= porder; ++k){
////
////				int nodei = Get_single_index(j, k, porder+1);
//				std::cout << j << " "  << " " << nodal::lagrange_l[i][j] << "\n";
////			}
//		}
////
//	}

	// terminate mpi
        int ierr = MPI_Finalize();	

}


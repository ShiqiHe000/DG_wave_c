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
#include "dg_verification.h"
#include <ctime>        // time()
#include <cstdlib>      // random numbe
#include "dg_cross_section.h"
//#include "dg_test.h"	//test


int main(int argc, char *argv[]){

	// initialize mpi
	Start_mpi(argc, argv);
	
	srand(time(NULL) + mpi::rank * mpi::num_proc);

	// prepare Hilbert curve
	Hilbert_numbering();
	
	// start parallel process
	Start_parallel();

	// start the game
	Driver_for_DG_approximation();

	// verification
	Get_error();
//	Write_error();	// write element-wise error to file

	// output cross section data
//	Solution_cross_section(0.55);


	// terminate mpi
        MPI_Finalize();	
	
}


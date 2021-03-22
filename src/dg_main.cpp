// Doxygen mainpage
/**
 * \author Shiqi He
 * \version 1.0
 * \mainpage  A Hash Table AMR Discontinuous Galerkin Spectral Element Wave Equation Solver with Space-Filling Curve Repartitioning Algorithm. 
 * \section intro_sec Introduction 
 * A hp-adaptive parallel Discontinuous Galerkin Spectral Element Solver
 * for an acoustic wave equation. The solver uses hash table data structure to 
 * arrange AMR data. A Space-filling curve (the HIlbert curve) is invoked for 
 * dynamic load-balancing. 
 * \section compile_sec Compilation
 * Please generate the Makefile by using cmake. 
 * \subsection prerequiste Prerequisite
 - g++ 7.5.0
 - cmake minimum requirement: 3.9 
 - OpenMPI: 4.0.2
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
#include "dg_total_element_num.h"
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

	// get total element num
//	Total_element_num(dg_time::t_total);

	// terminate mpi
        MPI_Finalize();	
	
}


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
#include "dg_read_mesh_2d.h"


int main(int argc, char *argv[]){
	
	int ierr;
	char** a = argv;	
	
	// Initialize mpi
	Start_mpi(&argc, a);
	
	// test read mesh	
    	Read_mesh_2d();

	// terminate mpi
        ierr = MPI_Finalize();	
}

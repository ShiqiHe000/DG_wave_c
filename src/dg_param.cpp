#include <iostream>
#include <mpi.h>
#include <string>
#include "dg_param.h"

//variable you must change--------------------------------------------

/// @brief mesh file information
/// @param fileplace The path of mesh file and mesh file name
namespace fileinfo{
	const std::string fileplace = "../gmsh_files/4_elements.msh";
}

/// @brief Domain size

/// @param exp_x element number in x direction(exponential order), i.e. 2^(exp_x)

/// @param exp_y element number in y direction(exponential order), i.e. 2^(exp_y)

/// @param gx_l domain left boundary (x direction, must >= 0) 
/// @param gx_r domain right boundary (x direction, must >= 0) 
/// @param gy_l domain left boundary (y direction, must >= 0) 
/// @param gy_r domain right boundary (y direction, must >= 0) 
/// @param nmin minimum polynomial degree in x direction
/// @param mmin minimum polynomial degree in y direction
/// @param nmax maximum polynomial degree in a direction
/// @param mmax maximum polynomial degree in y direction
namespace grid{
	const int exp_x = 1; 
	const int exp_y = 1; 
	
	const double gx_l = 0.0;
	const double gx_r = 1.0; 
	const double gy_l = 0.0;
	const double gy_r = 1.0; 

	const int nmin = 4;	// x direction
//	const int mmin = 4;	// y direction
	const int nmax = 10;
//	const int mmax = 10;

};
//---------------------------------------------------------------------


// variables you do not need to change----------------------------------

/// @brief mpi variables
/// @param rank process rank
/// @param num_proc total number of processor
namespace mpi{
	int rank;    
	int num_proc;    
}

//-----------------------------------------------------------------------

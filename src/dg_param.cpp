#include <string>
#include "dg_param.h"
#include <cmath>

//variable you must change--------------------------------------------

/// @brief mesh file information
/// @param fileplace The path of mesh file and mesh file name
namespace fileinfo{
//	const std::string fileplace = "../gmsh_files/sin_256.msh";
//	const std::string fileplace = "../gmsh_files/4_elem_small_domain.msh";
	const std::string fileplace = "../gmsh_files/128_128_mesh.msh";
}

/// @brief Domain size
/// @param exp_x element number in x direction(exponential order), i.e. 2^(exp_x)
/// @param exp_y element number in y direction(exponential order), i.e. 2^(exp_y)
/// @param gx_l domain left boundary (x direction, must >= 0) 
/// @param gx_r domain right boundary (x direction, must >= 0) 
/// @param gy_l domain left boundary (y direction, must >= 0) 
/// @param gy_r domain right boundary (y direction, must >= 0) 
/// @param nmin minimum polynomial degree in x and y direction
/// @param nmax maximum polynomial degree in x and y direction
/// @param hlevel_max maximum h-refinement level. 
namespace grid{
	const int exp_x = 7; 
	const int exp_y = 7; 
	
	const double gx_l = 0.0;
	const double gx_r = 16.0; 
	const double gy_l = 0.0;
	const double gy_r = 16.0; 

	const int nmin = 6;	
	const int nmax = 16;

	const int hlevel_max = 3;	
};
//---------------------------------------------------------------------


// variables you could change------------------------------------------

/// @brief
/// Time related variables
/// @param t_total total time integal
/// @param nt time step number
namespace dg_time{

	const double t_total = 1.0e-5 * 100;
//	const double t_total = 0.5;

	const int nt = 100;

};

/// @brief
/// Function related parameters
/// @param num_of_equation  number of equations
/// @param C wave speed
namespace dg_fun{

	const int num_of_equation = 3;

	const double C = 1.0;
//	const double C =  1.0 / (4.0 * std::atan(1.0));
};

/// @brief
/// refinement (hp) refinement swtich.
/// @param adapt refinement switch. 
/// @param refine_frequency time step interval of refinement. 
/// @param fit_point_num The number of points that are used to compute least square fit
/// @param tolerance_min the minimum discretization tolerance . If the estimated error exceeds this, refine. 
/// @param tolerance_max the maximum discretization tolerance. If the estimated error smaller than this, coarsen. 
/// @param load_balaning Repartitioning switch. 
namespace dg_refine{

	const bool adapt = true;

	const int refine_frequency = 20;	// every time step refine once

	const int fit_point_num = 4;

	const double tolerance_min = 1.0e-9;

	const double tolerance_max = 1.0e-14;

	const bool load_balancing = true;
};

//----------------------------------------------------------------------

/// @brief
/// variable to control the outputs
/// @param
namespace dg_io{

	const int output_frequency = 20;
};


// variables you do not need to change----------------------------------

/// @brief mpi variables
/// @param rank process rank
/// @param num_proc total number of processor
namespace mpi{
	int rank;    
	int num_proc;    
};

//-----------------------------------------------------------------------

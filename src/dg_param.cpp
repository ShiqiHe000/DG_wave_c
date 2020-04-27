#include <string>
#include "dg_param.h"
#include <cmath>

//variable you must change--------------------------------------------

/// @brief mesh file information
/// @param fileplace The path of mesh file and mesh file name
namespace fileinfo{
//	const std::string fileplace = "../gmsh_files/sin_256.msh";
//	const std::string fileplace = "../gmsh_files/4_elem_small_domain.msh";
	const std::string fileplace = "../gmsh_files/4_elements.msh";
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
	const int exp_x = 1; 
	const int exp_y = 1; 
	
	const double gx_l = 0.0;
	const double gx_r = 1.0; 
	const double gy_l = 0.0;
	const double gy_r = 1.0; 

	const int nmin = 8;	// x direction
	const int nmax = 14;

	const int hlevel_max = 3;	
};
//---------------------------------------------------------------------


// variables you could change------------------------------------------

/// @brief
/// Time related variables
/// @param t_total total time integal
/// @param nt time step number
namespace dg_time{

	const double t_total = 1.0e-5 * 2;
//	const double t_total = 0.5;

	const int nt = 2;

};

/// @brief
/// Function related parameters
/// @param num_of_equation  number of equations
/// @param C wave speed
namespace dg_fun{

	const int num_of_equation = 3;

//	const double C = 1.0;
	const double C =  1.0 / (4.0 * std::atan(1.0));
};

/// @brief
/// refinement (hp) refinement swtich.
/// @param adapt refinement switch. 
/// @param fit_point_num The number of points that are used to compute least square fit
/// @param tolerance the tolerance of global error. 
namespace dg_refine{

	const bool adapt = false;

	const int fit_point_num = 4;

	const double tolerance = 1.0e-9;

};

//----------------------------------------------------------------------

// variables you do not need to change----------------------------------

/// @brief mpi variables
/// @param rank process rank
/// @param num_proc total number of processor
namespace mpi{
	int rank;    
	int num_proc;    
};

//-----------------------------------------------------------------------

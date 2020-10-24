#ifndef DG_PARAM_H
#define DG_PARAM_H

#include <string>

// mesh file----------------------------------------
namespace fileinfo{
	extern const std::string fileplace;
	extern const std::string output_place;
	extern const std::string eff_filename;
	extern const std::string crosssection_filename;
	extern const std::string exact_error_filename;
	extern const std::string first_der_place;
	extern const std::string eigenvalues_place;
}
// -------------------------------------------------

// Grid size ---------------------------------------
namespace grid{
	extern const int exp_x;
	extern const int exp_y;

	extern const double gx_l;
	extern const double gx_r;
	extern const double gy_l;
	extern const double gy_r;
	
	extern const int nmin;
	extern const int nmax;

	extern const int hlevel_max;
};
//--------------------------------------------------

// time---------------------------------------------
namespace dg_time{

	extern const double t_total;
	extern const int nt;

};
//--------------------------------------------------


// dg function related parameters--------------------
namespace dg_fun{

	extern const int num_of_equation;
	
	extern const double C;

	extern const double pi;
};

//---------------------------------------------------

// refinemnt ----------------------------------------
namespace dg_refine{

	extern const bool adapt;

	extern const int refine_frequency;

	extern const int fit_point_num;
	
	extern const double tolerance_min;

	extern const double tolerance_max;

	extern const bool load_balancing;
};
//---------------------------------------------------

// output variables --------------------------------
namespace dg_io{

	extern const bool io;

	extern const int output_frequency;
};
//--------------------------------------------------

// mpi variables-------------------------------------
namespace mpi{
	extern int rank;
	extern int num_proc;
};
//---------------------------------------------------


#endif

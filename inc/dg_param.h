#ifndef DG_PARAM_H
#define DG_PARAM_H

#include <string>

// mesh file----------------------------------------
namespace fileinfo{
	extern const std::string fileplace;
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
};

//---------------------------------------------------

// refinemnt ----------------------------------------
namespace dg_refine{

	extern const bool adapt;

	extern const int fit_point_num;
};
//---------------------------------------------------


// mpi variables-------------------------------------
namespace mpi{
	extern int rank;
	extern int num_proc;
};
//---------------------------------------------------


#endif

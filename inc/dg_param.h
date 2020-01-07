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
};
//--------------------------------------------------

// time---------------------------------------------
namespace dg_time{

	extern const double t_total;
	extern const int nt;

};
//--------------------------------------------------

// mpi variables-------------------------------------
namespace mpi{
	extern int rank;
	extern int num_proc;
};
//---------------------------------------------------


#endif

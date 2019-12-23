#ifndef DG_PARAM_H
#define DG_PARAM_H

#include <iostream>
#include <string>

// mesh file----------------------------------------
namespace fileinfo{
	extern const std::string meshfile;
	extern const std::string fileplace;
}
// -------------------------------------------------


// mpi variables-------------------------------------
namespace mpi{
	extern int rank;
	extern int num_proc;
}
//---------------------------------------------------


#endif

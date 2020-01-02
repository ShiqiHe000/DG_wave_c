#include <mpi.h>
#include "dg_local_storage.h"

/// @brief Local data storage (on each process)
/// @param local_elem_num local element number, start with 1.
/// @param original_elem_num total element number before refine. 
/// @param elem_range the range of element stores on each processor.  Stored the last element numebr (start with 0, global index). 
/// @param x_local coordinate, start with 0.
/// @param y_local y coordinate, start with 0.
/// @param plevel_x polynomial level in x direction (start with 0)
/// @param plevel_y polynomial level in y direction (start with 0)
namespace local{
	
	int local_elem_num;

	int original_elem_num;

	int* elem_range;

	double* x_local = nullptr;
	double* y_local = nullptr;
	
	
	int* plevel_x = nullptr;
	int* plevel_y = nullptr;
};

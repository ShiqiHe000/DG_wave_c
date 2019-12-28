#include <mpi.h>

/// @brief Local data storage (on each process)
/// @param local_elem_num local element number, start with 1.
/// @param original_elem_num total element number before refine. 
/// @param elem_range the range of element stores on each processor. Stored the last element numebr (start with 0, global index). 
namespace local{
	
	int local_elem_num;

	int original_elem_num;

	int* elem_range;

};

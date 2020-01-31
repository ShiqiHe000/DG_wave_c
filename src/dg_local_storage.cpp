#include "dg_local_storage.h"
#include "dg_unit.h"	// head pointer
#include <unordered_map>
#include <vector>
#include "dg_boundary_table.h"

/// @brief Local data storage (on each process)
/// @param local_elem_num local element number, start with 1.
/// @param original_elem_num total element number before refine. 
/// @param elem_range the range of element stores on each processor.  Stored the last element numebr (start with 0, global index). 
/// @param rank_indicator store the last element Hilbert index on the finest level (useful in message exchange process).
/// @param x_local coordinate, start with 0.
/// @param y_local y coordinate, start with 0.
/// @param plevel_x polynomial level in x direction (start with 0)
/// @param plevel_y polynomial level in y direction (start with 0)
/// @param status element Hilbert status. 
/// @param head head pointer points to the first Unit. 
/// @param solution_int_l interior solution on the left boundary.
/// @param solution_int_r interior solution on the left boundary.
namespace local{
	
	int local_elem_num;

	int original_elem_num;	// no use?

	int* elem_range;
	int* rank_indicator; // no use?

	double* x_local = nullptr;
	double* y_local = nullptr;
	

	char* status;

	Unit* head = nullptr; // head ptr points to the first element

	std::unordered_map<int, Unit*> Hash_elem;	// hash table element

	std::vector<std::vector<double>> solution_int_l{};
	std::vector<std::vector<double>> solution_int_r{};

	// testing -----------------------------------------

//	std::vector<int> n_interface{};
//	std::vector<int> s_interface{};

	//--------------------------------------------------
};


/// @brief
/// The variables that will be used in h-refinement
namespace hrefinement{

	
	std::vector<table_elem> north;
	std::vector<table_elem> south;

	// accumulates boundaries info
	std::vector<accum_elem> south_accum;
	std::vector<accum_elem> north_accum;
};

/// @brief
/// Variables for results verification
namespace result{
	
	double* exact = nullptr;

	double* error = nullptr;

	double* L2_norm = nullptr;

};

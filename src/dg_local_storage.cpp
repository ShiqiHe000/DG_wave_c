#include "dg_local_storage.h"
#include "dg_unit.h"	// head pointer
#include <unordered_map>
#include <vector>
#include "dg_load_balancing.h"
#include "dg_boundary_table.h"

/// @brief Local data storage (on each process)
/// @param local_elem_num local element number, start with 1.
/// @param elem_range the range of element stores on each processor. Stored the last element numeber (start with 0, global index). 
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

	int* elem_range;

	double* x_local = nullptr;
	double* y_local = nullptr;
	

	char* status;

	Unit* head = nullptr; // head ptr points to the first element

	std::unordered_map<long long int, Unit*> Hash_elem;	// hash table element

};



/// @brief
/// The variables that will be used in h-refinement
/// @param north north mpi table. <target_rank, neighbour_keys>
/// @param south south mpi table.
/// @param west west mpi table.
/// @param east east mpi table.
/// @param neighbours_north Hash table to store all the possible neighbours in north direction. <key, all neighbours' key>
/// @param neighbours_south Hash table to store all the possible neighbours in south direction. <key, all neighbours' key>
/// @param neighbours_east Hash table to store all the possible neighbours in east direction. <key, all neighbours' key>
/// @param neighbours_west Hash table to store all the possible neighbours in west direction. <key, all neighbours' key>
namespace hrefinement{

	std::unordered_map<int, std::vector<mpi_table>> north;
	std::unordered_map<int, std::vector<mpi_table>> south;
	std::unordered_map<int, std::vector<mpi_table>> west;
	std::unordered_map<int, std::vector<mpi_table>> east;

	std::unordered_map<long long int, std::vector<long long int>> neighbours_north;
	std::unordered_map<long long int, std::vector<long long int>> neighbours_south;
	std::unordered_map<long long int, std::vector<long long int>> neighbours_east;
	std::unordered_map<long long int, std::vector<long long int>> neighbours_west;
};


/// @brief 
/// Variables for load balancing. 
namespace LB{
	
	std::vector<pmap> proc_mapping_table;

	double opt_bottleneck{};	// optimal bottleneck of the load balancing 

	int elem_accum{}; 	// elem prefix sum of former processor

	struct sending_envelope Send;	// record what to send
	
	Unit* my_rank_last = nullptr;	// pointer who points to the last element who stays in my rank. 
	Unit* my_rank_first = nullptr;	// pointer who points to the first element who stays in my rank. 

	bool first = true;	// first time evaluate the LB quality

	bool high_eff = false;	// repartitioning get 80% of the optimal efficiency
	
	double load_average{};	// ideal workload among processors (mean workload)
};

/// @brief
/// Variables for results verification
namespace result{
	
	std::unordered_map<int, std::vector<double>> exact;

	std::unordered_map<int, std::vector<double>> error;

	std::vector<double> L2_norm;

};

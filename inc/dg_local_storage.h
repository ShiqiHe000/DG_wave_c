#ifndef DG_LOCAL_STORAGE_H
#define DG_LOCAL_STORAGE_H

#include "dg_unit.h"
#include <unordered_map>
#include <vector>
#include "dg_load_struct.h"
#include "dg_boundary_table.h"

namespace local{
 
	extern int local_elem_num;

        extern int* elem_range;
	
	extern double* x_local; 
	extern double* y_local; 
	
	extern char* status;

	extern Unit* head;
	
	extern std::unordered_map<long long int, Unit*> Hash_elem;
	
};


namespace hrefinement{

	extern std::unordered_map<int, std::vector<mpi_table>> north;
	extern std::unordered_map<int, std::vector<mpi_table>> south;
	extern std::unordered_map<int, std::vector<mpi_table>> west;
	extern std::unordered_map<int, std::vector<mpi_table>> east;

	
	extern std::unordered_map<long long int, std::vector<long long int>> neighbours_north;
	extern std::unordered_map<long long int, std::vector<long long int>> neighbours_south;
	extern std::unordered_map<long long int, std::vector<long long int>> neighbours_east;
	extern std::unordered_map<long long int, std::vector<long long int>> neighbours_west;
};

namespace LB{

	extern std::vector<pmap> proc_mapping_table;
	
	extern double opt_bottleneck;

	extern int elem_accum;

	extern struct sending_envelope Send;
	
	extern Unit* my_rank_last;
	extern Unit* my_rank_first;
	
	extern bool first;

	extern bool high_eff;
	
	extern double load_average;
};

namespace result{

	extern std::unordered_map<int, std::vector<double>> exact;

	extern std::unordered_map<int, std::vector<double>> error;

	extern std::vector<double> L2_norm;
};

#endif

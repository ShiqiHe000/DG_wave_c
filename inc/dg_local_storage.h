#ifndef DG_LOCAL_STORAGE_H
#define DG_LOCAL_STORAGE_H

#include "dg_unit.h"
#include <unordered_map>
#include <vector>
#include "dg_boundary_table.h" 
#include "dg_load_struct.h"

namespace local{
 
	extern int local_elem_num;
 
        extern int original_elem_num;

        extern int* elem_range;
	extern int* rank_indicator;
	
	extern double* x_local; 
	extern double* y_local; 
	
	extern char* status;

	extern Unit* head;
	
	extern std::unordered_map<int, Unit*> Hash_elem;

	extern std::vector<std::vector<double>> solution_int_l;
	extern std::vector<std::vector<double>> solution_int_r;
	
};


namespace hrefinement{

	extern std::unordered_map<int, std::vector<mpi_table>> north;
	extern std::unordered_map<int, std::vector<mpi_table>> south;
	extern std::unordered_map<int, std::vector<mpi_table>> west;
	extern std::unordered_map<int, std::vector<mpi_table>> east;

	
	extern std::unordered_map<int, std::vector<int>> neighbours_north;
	extern std::unordered_map<int, std::vector<int>> neighbours_south;
	extern std::unordered_map<int, std::vector<int>> neighbours_east;
	extern std::unordered_map<int, std::vector<int>> neighbours_west;
};

namespace LB{

	extern std::vector<pmap> proc_mapping_table;

	extern int elem_accum;

	extern struct sending_envelope Send;
	
	extern Unit* end;
	extern Unit* my_rank_last;
};

namespace result{

	extern double* exact;
	
	extern double* error;

	extern double* L2_norm;

};

#endif

#ifndef DG_LOCAL_STORAGE_H
#define DG_LOCAL_STORAGE_H

#include "dg_unit.h"
#include <unordered_map>
#include <vector>

namespace local{
 
	extern int local_elem_num;
 
        extern int original_elem_num;
 
        extern int* elem_range;
	
	extern double* x_local; 
	extern double* y_local; 
	
	extern char* status;

	extern Unit* head;
	
	extern std::unordered_map<int, Unit*> Hash_elem;

	extern std::vector<std::vector<double>> solution_int_l;
	extern std::vector<std::vector<double>> solution_int_r;
	
	extern std::vector<std::vector<double>> ghost;
};

namespace result{

	extern double* exact;
	
	extern double* error;

	extern double* L2_norm;

};

#endif

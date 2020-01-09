#ifndef DG_LOCAL_STORAGE_H
#define DG_LOCAL_STORAGE_H

#include "dg_unit.h"

namespace local{
 
	extern int local_elem_num;
 
        extern int original_elem_num;
 
        extern int* elem_range;
	
	extern double* x_local; 
	extern double* y_local; 
	
	extern int* plevel_x;
	extern int* plevel_y; 

	extern char* status;

	extern Unit* head;
};


#endif

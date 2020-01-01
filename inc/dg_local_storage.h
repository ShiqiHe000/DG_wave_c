#ifndef DG_LOCAL_STORAGE_H
#define DG_LOCAL_STORAGE_H

namespace local{
 
	extern int local_elem_num;
 
        extern int original_elem_num;
 
        extern int* elem_range;
	
	extern double* x_local; 
	extern double* y_local; 
	
	extern double* pl_p[];
	extern double* pl_w[];
};


#endif

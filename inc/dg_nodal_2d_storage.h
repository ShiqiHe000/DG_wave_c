#ifndef DG_NODAL_2D_STORAGE_H
#define DG_NODAL_2D_STORAGE_H

namespace SortMesh{

	extern double* elem_x_position;
	extern double* elem_y_position;
	
	extern int num_of_element_x;
	extern int num_of_element_y;
	extern int num_of_element;
	
	extern int* dual_coord;
	
	extern double* x_hilbert;
	extern double* y_hilbert;
}

namespace nodal{

	extern double* gl_p[];
	extern double* gl_w[];
	
	extern double* first_der[];
}

#endif

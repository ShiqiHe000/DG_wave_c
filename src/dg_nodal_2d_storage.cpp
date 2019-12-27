#include <mpi.h>

/// @brief vairables used in sort mesh file
/// @param elem_x_position kth element x coordinates in design sequence
/// @param elem_y_position kth element y coordinates in design sequence
/// @param num_of_element_x number of element in x direction (start with 1)
/// @param num_of_element_y number of element in y direction (start with 1)
/// @param num_of_element Total number of element (originally)
/// @param dual_coord dual graph coordinates
/// @param x_hilbert x coordinate of element k on hilbert curve
/// @param y_hilbert y coordinate of element k on hilbert curve
namespace SortMesh{
	double* elem_x_position; 
	double* elem_y_position; 

	int num_of_element_x;	
	int num_of_element_y;	
	int num_of_element;
	
	int* dual_coord; 
	
	double* x_hilbert;
	double* y_hilbert;
}


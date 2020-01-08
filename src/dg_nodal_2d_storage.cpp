#include <mpi.h>
#include "dg_nodal_2d_storage.h"

/// @brief vairables used in sort mesh file
/// @param elem_x_position kth element x coordinates in design sequence
/// @param elem_y_position kth element y coordinates in design sequence
/// @param num_of_element_x number of element in x direction (start with 1)
/// @param num_of_element_y number of element in y direction (start with 1)
/// @param num_of_element Total number of element (originally)
/// @param dual_coord dual graph coordinates
/// @param x_hilbert x coordinate of element k on hilbert curve
/// @param y_hilbert y coordinate of element k on hilbert curve
/// @param status element status
namespace SortMesh{
	double* elem_x_position = nullptr; 
	double* elem_y_position = nullptr; 

	int num_of_element_x;	
	int num_of_element_y;	
	int num_of_element;
	
	int* dual_coord; 
	
	double* x_hilbert = nullptr;
	double* y_hilbert = nullptr;

	char* status = nullptr;
}

/// @brief 
/// 2d nodal data storage
/// @param gl_p arrau of pointer, each pointer points to an array of GL points of corresponding poly order.
/// @param gl_w arrau of pointer, each pointer points to an array of GL weights of corresponding poly order.
/// @param first_der first derivative matrix.
/// @param lagrange_l Lagrange interpolating polynomial on the left boundary
/// @param lagrange_r Lagrange interpolating polynomial on the right boundary
namespace nodal{

	double** gl_p{};
	double** gl_w{};

	double** first_der{};
		
	double** lagrange_l{};
	double** lagrange_r{};
}

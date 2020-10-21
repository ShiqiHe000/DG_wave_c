#include "dg_nodal_2d_storage.h"
#include <vector>
#include <unordered_map>

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

	std::vector<char> status; 
}

/// @brief 
/// 2d nodal data storage
/// @param gl_points Hash table for Gauss points. Key: polynoial order, value: vector of guass points.
/// @param gl_weights Hash table for Gauss weights. Key: polynoliam order, value: Gauss weigths. 
/// @param first_der First derivative matrix. Key: poynomial order. 
/// @param mfirst_der Modified first derivative matrix. Key: poynomial order. 
/// @param lagrange_l Lagrange interpolating polynomial on the left boundary. Key: poynomial order.
/// @param lagrange_r Lagrange interpolating polynomial on the right boundary. Key: poynomial order.
/// @param legendre Legendre polynomial value at GL points. Key: polynomial order. 
/// @param gll_points Hash table for Gauss Lobbato points. Key: polynoial order, value: vector of guass points.
/// @param gll_weights Hash table for Gauss Lobbato weights. Key: polynoliam order, value: Gauss weigths. 
namespace nodal{

	std::unordered_map<int, std::vector<double>> gl_points;
	std::unordered_map<int, std::vector<double>> gl_weights;

	std::unordered_map<int, std::vector<double>> first_der;
	std::unordered_map<int, std::vector<double>> mfirst_der;

	std::unordered_map<int, std::vector<double>> lagrange_l;
	std::unordered_map<int, std::vector<double>> lagrange_r;

	std::unordered_map<int, std::vector<double>> legendre;	// useless

	std::unordered_map<int, std::vector<double>> gll_points;
	std::unordered_map<int, std::vector<double>> gll_weights;
}

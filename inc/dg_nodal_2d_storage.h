#ifndef DG_NODAL_2D_STORAGE_H
#define DG_NODAL_2D_STORAGE_H

#include <vector>
#include <unordered_map>

namespace SortMesh{

	extern double* elem_x_position;
	extern double* elem_y_position;
	
	extern int num_of_element_x;
	extern int num_of_element_y;
	extern int num_of_element;
	
	extern int* dual_coord;
	
	extern double* x_hilbert;
	extern double* y_hilbert;

	extern std::vector<char> status;
}

namespace nodal{

//	extern double** gl_p;
//	extern double** gl_w;

	extern std::unordered_map<int, std::vector<double>> gl_points;
	extern std::unordered_map<int, std::vector<double>> gl_weights;

	extern std::unordered_map<int, std::vector<double>> first_der;

	extern std::unordered_map<int, std::vector<double>> lagrange_l;
	extern std::unordered_map<int, std::vector<double>> lagrange_r;

//	extern double** first_der;
//
//	extern double** lagrange_l;
//	extern double** lagrange_r;

}

#endif

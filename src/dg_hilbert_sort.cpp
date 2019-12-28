#include <iostream>
#include <mpi.h>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_hilbert_curve.h"
#include "dg_single_index.h"

void Hilbert_sort_2d(){

	// economic storage (only store the 2 diagonal nodes)
	SortMesh::x_hilbert = new double[2 * SortMesh::num_of_element]{};
	SortMesh::y_hilbert = new double[2 * SortMesh::num_of_element]{};

	for(int k = 0; k < SortMesh::num_of_element; ++k){
	
		int i = SortMesh::dual_coord[Get_single_index(k, 0, 2)];
		int j = SortMesh::dual_coord[Get_single_index(k, 1, 2)];

	
		int d = xy2d(grid::exp_x, j, i);
		
		int node1 = Get_single_index(d, 0, 2);
		int node2 = Get_single_index(d, 1, 2);
		SortMesh::x_hilbert[node1] = SortMesh::elem_x_position[Get_single_index(k, 0, 4)]; 
		SortMesh::x_hilbert[node2] = SortMesh::elem_x_position[Get_single_index(k, 2, 4)]; 
		SortMesh::y_hilbert[node1] = SortMesh::elem_y_position[Get_single_index(k, 0, 4)]; 
		SortMesh::y_hilbert[node2] = SortMesh::elem_y_position[Get_single_index(k, 2, 4)]; 

	}
	
	// free memory--------------------------------------
	delete[] SortMesh::elem_x_position;
	delete[] SortMesh::elem_y_position;
	delete[] SortMesh::dual_coord;

	SortMesh::elem_x_position = nullptr;
	SortMesh::elem_y_position = nullptr;
	SortMesh::dual_coord = nullptr;
	//--------------------------------------------------

	std::cout<< "-----------------------------------------" << "\n";
	std::cout<< "Finished read mesh file and sorting." << "\n";
	std::cout<< "-----------------------------------------" << "\n";
}

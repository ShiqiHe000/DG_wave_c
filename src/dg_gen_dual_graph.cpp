#include <mpi.h>
#include <cmath>	// pow()
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"

void Gen_dual_graph_2d(){
	
	// compute the domain size----------------------------------------------------------
	SortMesh::num_of_element_x = std::pow(2, grid::exp_x);		
	SortMesh::num_of_element_y = std::pow(2, grid::exp_y);		
	SortMesh::num_of_element = SortMesh::num_of_element_x * SortMesh::num_of_element_y;
	//---------------------------------------------------------------------------------
	
	// allocate---------------------------------------------------------
	SortMesh::dual_coord = new int[SortMesh::num_of_element * 2]{};
	//-------------------------------------------------------------------

	// element size (start with uniform size)-------------------------
	double delta_x = (grid::gx_r - grid::gx_l) / (double)SortMesh::num_of_element_x;
	double delta_y = (grid::gy_r - grid::gy_l) / (double)SortMesh::num_of_element_y;
	//----------------------------------------------------------------

	for(int i = 0; i < SortMesh::num_of_element; ++i){
		

	}

	

	

}

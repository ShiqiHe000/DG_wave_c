#include <mpi.h>
#include <cmath>	// pow()
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_get_dual_coord.h"
#include "dg_single_index.h"
#include "dg_gen_dual_graph.h"

/// @brief
/// Generate the dual graph coordinates of the mesh file. n(i, j).
/// Coordinate expansion: 0 <= i <= EXP_X - 1, 0 <= j <= EXP_Y - 1
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

	for(int k = 0; k < SortMesh::num_of_element; ++k){
		int inode = Get_single_index(k, 0, 4);
		int inode2 = Get_single_index(k, 0, 2);
		int* temp = &SortMesh::dual_coord[inode2];
		Get_dual_coord_2d(SortMesh::elem_x_position[inode], SortMesh::elem_y_position[inode], 
					delta_x, delta_y, temp);		

	}

	
	

}

#include <mpi.h>
#include <cmath>	// std::pow
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include "dg_distribute_elements.h"
#include "dg_start_parallel.h"

void Start_parallel(){
	
	// Distribute the elements evenly between processors
	Distribute_elem();

	SortMesh::num_of_element_x = std::pow(2, grid::exp_x); 
	SortMesh::num_of_element_y = std::pow(2, grid::exp_y); 



}

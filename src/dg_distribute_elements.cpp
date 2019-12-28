#include <mpi.h>
#include <algorithm>	// std::fill_n
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include "dg_local_storage.h"
#include "dg_distribute_elements.h"
#include <iostream> 	// just for test

/// @brief
/// After reading the mesh file, we spread the elements between processors 
/// equally. We support element number < processors number, i.e., a 
/// processor could be assigned with zero elements.
void Distribute_elem(){
	
	// Integer array (of length group size) specifying the number of elements to send to each processor.
	int sendcouts[mpi::num_proc]{};
	
	// Integer array (of length group size). 
	//Entry i specifies the displacement (relative to sendbuf) from which to take the outgoing data to process i. 
	int displs[mpi::num_proc]{};
	
	// allocate one more unit for the offset
	local::elem_range = new int[mpi::num_proc + 1]{};
	local::elem_range[0] = -1;
	
	// distribute elements
	int local_elem_number[mpi::num_proc]{}; 
	
	// average load
	int average;

	// last processor's load
	int last;

	// process 1 compute the average load
	if(mpi::rank == 0){
		
		local::original_elem_num = SortMesh::num_of_element;

		// if element number < processor number
		if(SortMesh::num_of_element < mpi::num_proc){
			
			for(int ii = 0; ii < SortMesh::num_of_element; ++ii){
				local_elem_number[ii] = 1;
				sendcouts[ii] = 1;
	
			}
				local::elem_range[1] = 0;
		}
		else{	// element number >= processor number
			average = SortMesh::num_of_element / mpi::num_proc;
			last = SortMesh::num_of_element - (mpi::num_proc -1) * average;

			local::elem_range[1] = average - 1;	// element numbering start with 0
			
			std::fill_n(local_elem_number, mpi::num_proc - 1, average);
			
			local_elem_number[mpi::num_proc - 1] = last;

			std::fill_n(sendcouts, mpi::num_proc - 1, average);

			sendcouts[mpi::num_proc -1 ] = last; 

		}

	}
}

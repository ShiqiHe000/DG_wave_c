#include <mpi.h>
#include <algorithm>	// std::fill_n
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include "dg_local_storage.h"
#include "dg_distribute_elements.h"
#include <cmath>	// pow
#include <vector>
#include <iostream>	//test

/// @brief
/// After reading the mesh file, we spread the elements between processors 
/// equally. We support element number < processors number, i.e., a 
/// processor could be assigned with zero elements.
void Distribute_elem(){
	
	// Integer array (of length group size) specifying the number of elements to send to each processor.
	int sendcouts[mpi::num_proc]{};	// each element 2 x/y coordinates (diagonal). 
	int sendcouts_status[mpi::num_proc]{};	// each element 1 status
	
	// Integer array (of length group size). 
	//Entry i specifies the displacement (relative to sendbuf) from which to take the outgoing data to process i. 
	int displs[mpi::num_proc]{};
	int displs_status[mpi::num_proc]{};
	
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
		
		// if element number < processor number
		if(SortMesh::num_of_element < mpi::num_proc){
			
			for(int ii = 0; ii < SortMesh::num_of_element; ++ii){
				local_elem_number[ii] = 1;
				sendcouts[ii] = 1 * 2;	// for each element we record two diagnoal points
				sendcouts_status[ii] = 1;	
			}
				local::elem_range[1] = 0;
		}
		else{	// element number >= processor number
//std::cout<< "num_elem "<< SortMesh::num_of_element<< "\n";
			average = SortMesh::num_of_element / mpi::num_proc;
			last = SortMesh::num_of_element - (mpi::num_proc -1) * average;

			local::elem_range[1] = average - 1;	// element numbering start with 0
			
			std::fill_n(local_elem_number, mpi::num_proc - 1, average);		
	
			local_elem_number[mpi::num_proc - 1] = last;

			std::fill_n(sendcouts, mpi::num_proc - 1, average * 2);
			std::fill_n(sendcouts_status, mpi::num_proc - 1, average);

			sendcouts[mpi::num_proc -1 ] = last * 2; 
			sendcouts_status[mpi::num_proc -1 ] = last; 

		}

		for(int i = 1; i < mpi::num_proc; ++i ){
			
			displs[i] = local_elem_number[i-1] * 2 + displs[i-1];
			displs_status[i] = local_elem_number[i-1] + displs_status[i-1];
			local::elem_range[i+1] = local::elem_range[i] + local_elem_number[i];
		}

	}
	
	// scatter local element number
	MPI_Scatter(&local_elem_number[0], 1, MPI_INT, &local::local_elem_num, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	
	// allocate local storage
	local::x_local = new double[2 * local::local_elem_num ];
	local::y_local = new double[2 * local::local_elem_num ];
	local::status = new char[local::local_elem_num];

	// scatter data
	MPI_Scatterv(SortMesh::x_hilbert, sendcouts, displs, MPI_DOUBLE, local::x_local, local::local_elem_num * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(SortMesh::y_hilbert, sendcouts, displs, MPI_DOUBLE, local::y_local, local::local_elem_num * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(&SortMesh::status[0], sendcouts_status, displs_status, MPI_CHAR, local::status, local::local_elem_num, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(local::elem_range, mpi::num_proc + 1, MPI_INT, 0, MPI_COMM_WORLD);

	
	// deallocate
	if(mpi::rank == 0){
		delete[] SortMesh::x_hilbert;
		delete[] SortMesh::y_hilbert;
		//delete[] SortMesh::status;
		SortMesh::status.clear();

		SortMesh::x_hilbert = nullptr;
		SortMesh::y_hilbert = nullptr;
		//SortMesh::status = nullptr;
	}
}

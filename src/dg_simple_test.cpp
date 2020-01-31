#include "dg_simple_test.h"
#include "dg_unit.h"
#include <vector>
#include "dg_local_storage.h"
#include "dg_elem_length.h"
#include <mpi.h>
#include "dg_cantor_pairing.h"
#include "dg_boundary_table.h"
#include  <iostream>	// test
#include "dg_param.h"	//test


// forward declaration----------------------------------------------------------
void Boundary_condition(double& s_interface, int nt);
//------------------------------------------------------------------------------


/// @parma tn nth-delta_t.
void Simple_test(int tn){

	// tarverse the hash table
	Unit* temp = local::head;
	
	// total number of info to send
	int mpi_num{};
	for(auto& v : hrefinement::north_accum){
		mpi_num += v.sum;

	}


	MPI_Request request[mpi_num];
	MPI_Status status[mpi_num];
	int i{};

	// form interface + send info on the mpi boundary
	for(int k = 0; k < local::local_elem_num; ++k){

		temp -> n_interface = (double)(temp -> var);

		// send if on the mpi boudary
		for(auto& v : temp -> facen[1]){

			if(v.face_type == 'M'){
				int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				//tag=sender's key
				MPI_Isend(&(temp -> n_interface), 1, MPI_DOUBLE, v.rank, local_key, MPI_COMM_WORLD, &request[i]); 
//if(mpi::rank == 0){
//	std::cout<< "local_key " << local_key << " n_interface: " << temp -> n_interface << "\n";
//
//}

				++i;
			}
		}

		temp = temp -> next;
	}

	int length_tol = Elem_length(0);

	temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){
		
		if(temp -> facen[0].front().face_type == 'B'){	// if on the physical boundary
				
			Boundary_condition(temp -> s_interface, tn);

		}
		else{	// local/MPI neighbour

		 	for(auto& v : temp -> facen[0]){

				if(v.face_type == 'L'){	// local
					
					if(v.hlevel <= (temp -> index[2])){	// if neighbour is larger than current element, or equal size

						temp -> s_interface = local::Hash_elem[v.key] -> n_interface;

					}
					else{	// neighbour is smaller
					
						int n_length = Elem_length(v.hlevel);
						temp -> s_interface += (local::Hash_elem[v.key] -> n_interface * 
										(double)(n_length) / length_tol);
					}
				}
				else{	// neighbour is stored remotely

					MPI_Status status;
				
					int recv_info;

					MPI_Recv(&recv_info, 1, MPI_DOUBLE, v.rank, v.key, MPI_COMM_WORLD, &status); // tag = n_key
						
					if(v.hlevel <= (temp -> index[2])){	// if remote element is larger

						temp -> s_interface += recv_info;						

					}
					else{	// remote element is smaller

						
						int n_length = Elem_length(v.hlevel);
						temp -> s_interface += recv_info * (double)(n_length) / length_tol;

					}	

				}

			}


		}


		temp = temp -> next;

	}
	
	// confirm all the messages being received
	MPI_Waitall(mpi_num, request, &status[0]);

	// updates var
	temp = local::head;
	for(int k = 0; k < local::local_elem_num; ++k){

		temp -> var = (int)(temp -> s_interface);
		

		temp = temp -> next;
	}


}

/// @breif
/// Boundary condition. Imposed on the south interfaces.
/// @param s_interface south interfaces vector.
/// @param nt nth-delta_t
void Boundary_condition(double& s_interface, int nt){

	if((nt % 2) == 0){	// even time step

		s_interface = 1;

	}
	else{	// odd time step
		s_interface = -1;

	}

}

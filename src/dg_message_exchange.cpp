#include "dg_message_exchange.h"
#include "dg_local_storage.h"
#include "dg_unit.h"

/// @brief 
/// exchange information on the mpi boundary
void Exchange_solution_x(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		

		temp = temp->next;
	}
	

}

/// @brief
/// Inside one element use non-blocking send and blocking recv
/// to update MPI interfaces information. x direction. 
void Non_blocking_exchange_interface_x(Unit*& temp){
	
	//now only support conforming interface
	if((temp->faces[0]) == 1){	// face 1 recv



		int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
    int source, int tag, MPI_Comm comm, MPI_Status *status)

	}


}

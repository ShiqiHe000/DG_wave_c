#include <mpi.h>
#include "dg_param.h"
#include "dg_read_mesh_2d.h"


void Hilber_numbering(){
	
	if(mpi::rank == 0){
	
		Read_mesh_2d();		


	}


}

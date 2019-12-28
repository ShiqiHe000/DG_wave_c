#include <mpi.h>
#include "dg_param.h"
#include "dg_read_mesh_2d.h"
#include "dg_gen_dual_graph.h"
#include "dg_hilbert_sort.h"

void Hilber_numbering(){
	
	if(mpi::rank == 0){
		
		// read .msh file
		Read_mesh_2d();		
		
		// generate dual graph
		Gen_dual_graph_2d();
		
		// renumbering elements by hilbert curve
		Hilbert_sort_2d();

	}


}

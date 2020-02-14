#include "dg_adapt.h"
#include "dg_hrefinement.h"
#include "dg_local_storage.h"
#include "dg_mpi_table_construct.h"

void Adapt(){

	h_refinement();

	Construct_mpi_table(hrefinement::north, hrefinement::south);
	Update_mpi_boundaries(hrefinement::north, hrefinement::south);	

}

/// @brief
/// Remove all the elements inside the tables (MPI tables and accumlation table)
void Clear_tables(){

	hrefinement::south_accum.clear();
	hrefinement::north_accum.clear();
	hrefinement::south.clear();
	hrefinement::north.clear();


}

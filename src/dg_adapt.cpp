#include "dg_adapt.h"
#include "dg_hrefinement.h"
#include "dg_local_storage.h"
#include "dg_mpi_table_construct.h"

void Adapt(){

	h_refinement();
	
	// x direction
	Construct_mpi_table(hrefinement::north, 1, hrefinement::south, 0);
	Update_mpi_boundaries(hrefinement::north, 1, hrefinement::north_accum, hrefinement::south, 0, hrefinement::south_accum);	

	// y direction
	Construct_mpi_table(hrefinement::west, 2, hrefinement::east, 3);
	Update_mpi_boundaries(hrefinement::west, 2, hrefinement::west_accum, hrefinement::east, 3, hrefinement::east_accum);	
}

/// @brief
/// Remove all the elements inside the tables (MPI tables and accumlation table)
void Clear_tables(){

	hrefinement::south_accum.clear();
	hrefinement::north_accum.clear();
	hrefinement::south.clear();
	hrefinement::north.clear();

	hrefinement::west_accum.clear();
	hrefinement::east_accum.clear();
	hrefinement::west.clear();
	hrefinement::east.clear();

}

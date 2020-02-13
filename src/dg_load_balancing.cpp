#include "dg_load_balancing.h"
#include "dg_reallocate_elem.h"
#include "dg_proc_mapping.h"

/// @brief
/// Whole procedure of load balancing   
void Load_balancing(){

	Build_mapping_table();

	Update_mpi_boundary();

	Reallocate_elem();

}

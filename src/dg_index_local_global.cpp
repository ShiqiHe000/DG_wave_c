#include "dg_index_local_global.h"
#include "dg_local_storage.h"

/// @brief
/// Input rank and element local index, output global index.
/// @param rank process rank
/// @param l_index element local index
int Index_local_to_global(int rank, int l_index){
	
	return (l_index + local::elem_range[rank] + 1);

}

/// @brief
/// Input rank and element global index, output the local index.
/// @param rank process rank
/// @param g_index element global index
int Index_global_to_local(int rank, int g_index){

	return (g_index - (local::elem_range[rank] + 1));

}

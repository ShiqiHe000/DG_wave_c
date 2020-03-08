#ifndef DG_BOUNDARY_TABLE_H
#define DG_BOUNDARY_TABLE_H

/// @brief
/// Structure for MPI table. 
/// @param local_key Current element's key.
/// @param mpi_length The length of current element that is exposed to target_rank.
/// @param owners_rank The rank of the owner of this element after repartitioning. 
struct mpi_table{

	int local_key;

	int mpi_length;

	int owners_rank;
};



#endif



#ifndef DG_BOUNDARY_TABLE_H
#define DG_BOUNDARY_TABLE_H

/// @brief
/// MPI boundary information table.
/// @param local_key the key of current elememt.
/// @param target_rank neighbour's rank.
/// @param coord element coordinate of related direction. (x:j, y:i).
/// @param hlevel current element's hlevel.
/// @param mpi_length The portion of the element interface length that exposes to the mpi boundary. 
/// @param owners_rank The future (after repartitioning) rank of current element.
struct table_elem{

	int local_key;

	int target_rank;

	double coord;

	int hlevel;

	int mpi_length;
	
	int owners_rank;
};


/// @brief
/// Accumulation table. Store the neighbour ranks and the number of elements that is facing to the corresponding rank.
/// @param rank neighbour rank.
/// @param sum element number in current rank that is facing to the neighbour rank.
struct accum_elem{

	int rank;
	
	int sum{};

};



#endif



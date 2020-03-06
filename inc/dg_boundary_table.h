#ifndef DG_BOUNDARY_TABLE_H
#define DG_BOUNDARY_TABLE_H

/// @brief
/// Structure for MPI table. 
/// @param local_key Current element's key.
/// @param target_rank Neighbour's rank.
/// @param hlevel Current element's h-refinemnt level.
/// @param mpi_length The length of current element that is exposed to target_rank.
/// @param owners_rank The rank of the owner of this element after repartitioning. 
struct mpi_table{

	int local_key;

	int target_rank;

	int hlevel;

	int mpi_length;

	int owners_rank;
};

/// @brief
/// One element could have maximum two remote nieghbours in one direction (x or y).
/// Store the possible neighbours' key the corresponding vector. Left == face 0 or face 2.
/// Right == face 1 for face 3.
//struct two_sides{
//
//	std::vector<int> left;
//	
//	std::vector<int> right;
//};


#endif



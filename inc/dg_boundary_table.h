#ifndef DG_BOUNDARY_TABLE_H
#define DG_BOUNDARY_TABLE_H

/// @brief
/// MPI boundary information table.
/// @param local_key the key of current elememt.
/// @param target_rank neighbour's rank.
/// @param coord element coordinate of related direction. (x:j, y:i)
struct table_elem{

	int local_key;

	int target_rank;

	int coord;

	int hlevel;
};

/// @brief
/// Sorting function. element in coord ascending sequence.
bool compare_coord(table_elem left, table_elem right) {return left.coord < right.coord;}

struct accum_elem{

	int rank;
	
	int sum{};

};

#endif



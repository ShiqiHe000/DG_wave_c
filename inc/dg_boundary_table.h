#ifndef DG_BOUNDARY_TABLE_H
#define DG_BOUNDARY_TABLE_H

/// @brief
/// MPI boundary information table.
/// @param local_key the key of current elememt.
/// @param target_rank neighbour's rank.
/// @param coord element coordinate of related direction. (x:j, y:i).
/// @param hlevel current element's hlevel.
struct table_elem{

	int local_key;

	int target_rank;

	double coord;

	int hlevel;
};

bool compare_coord(table_elem left, table_elem right); 

struct accum_elem{

	int rank;
	
	int sum{};

};

///// @brief
///// compute the element side length. length = pow(2, hlevel_max - level)
//int Elem_length(int level){
//
//	return ((int)(std::pow(2, grid::hlevel_max - level) + 0.5));
//}

#endif



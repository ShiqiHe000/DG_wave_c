#include <mpi.h>
#include "dg_poly_level_and_order.h"

/// @brief
/// Transform polynomial order to polynomial level. Level start with 0.
/// @param nmin minimum polynomial order
/// @param plevel current element polynomial level
int Poly_level_to_order(int nmin, int plevel){

	return (nmin + 2 * plevel);

}

/// @brief 
/// Transform polynomial order to level
/// @param nmin minimum polynomial order
/// @param porder cuurent element polynomial order
int Poly_order_to_level(int nmin, int porder){

	return (porder - nmin) / 2;

}




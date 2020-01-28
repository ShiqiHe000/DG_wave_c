#include "dg_elem_length"
#include "dg_param.h"
#include <cmath>

/// @brief
/// compute the element side length. length = pow(2, hlevel_max - level)
int Elem_length(int level){

	 return ((int)(std::pow(2, grid::hlevel_max - level) + 0.5)); 

}

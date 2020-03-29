#include "dg_affine_map.h"
/// @brief
/// We map each element to the reference interval [-1, 1]. which will 
/// serve as our refernece element.
/// @param xi reference location.
/// @param xk1 the location of the left bound of this element.
/// @param delta_x element size
double Affine_mapping(double xi, double xk1, double delta_x){

	return (xk1 + (xi + 1.0) * delta_x / 2.0);

}

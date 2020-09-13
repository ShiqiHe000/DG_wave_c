#include "dg_affine_map.h"
/// @brief
/// We map each element to the reference interval [-1, 1]. which will 
/// serve as our refernece element. 
/// This function map element reference location back to phyical location. 
/// Input: the reference location. Output: the coordinate back to original coordinate. 
/// @param xi reference location.
/// @param xk1 the location of the left bound of this element.
/// @param delta_x element size
double Affine_mapping(double xi, double xk1, double delta_x){

	return (xk1 + (xi + 1.0) * delta_x / 2.0);

}

/// @brief
/// Map the phyical location to reference space. Output: coordinate on reference space. 
/// @param x physical location. 
/// @param x_l element left boundary coordinate. 
/// @param del_x element size. 
double Map_to_reference(double x, double x_l, double del_x){

	return ((x - x_l) / del_x * 2.0 - 1.0);

}

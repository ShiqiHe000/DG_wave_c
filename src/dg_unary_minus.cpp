#include <vector>
#include "dg_unary_minus.h"
#include <cassert>
#include <algorithm>

/// @brief
/// Vector b = - vector a (same size)
/// @param a the vector to apply - operator.
/// @param b vector b gets the value.
void Unary_minus(std::vector<double>& a, std::vector<double>& b){

	assert(a.size() == b.size() && "The size of the two vector should be equal. ");

	std::transform(a.begin(), a.end(), b.begin(), [](double x){return -x;});
}

#include <vector>
#include "dg_vector_operation.h"
#include <cassert>
#include <algorithm>
#include <functional>

/// @brief
/// Vector b = - vector a (same size)
/// @param a the vector to apply - operator.
/// @param b vector b gets the value.
void Unary_minus(std::vector<double>& a, std::vector<double>& b){

	assert(a.size() == b.size() && "The size of the two vector should be equal. ");

	std::transform(a.begin(), a.end(), b.begin(), [](double x){return -x;});
}


/// @brief
/// Performs the - operator between two vectors (b = b - a). The size of the two vector should equal.
/// @param a the vector of be substructed.
/// @param b the vector to be substructed from.
void Vector_minus(std::vector<double>& a, std::vector<double>& b){

	assert(a.size() == b.size() && "The size of the two vector should be equal. ");

	std::transform(b.begin(), b.end(), a.begin(), b.begin(), std::minus<double>());
}

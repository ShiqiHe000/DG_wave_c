#include "dg_unit.h"
#include "dg_param.h"
#include <vector>
#include <algorithm>


/// @brief
/// constructor (default).
/// default: level = 0, n = min poly order, m = max poly order
/// index = {0, 0, 0}, mpi_f = false, x/ycoords = {0.0, 0.0}
/// solution = nullptr. 
Unit::Unit() : n(grid::nmin), m(grid::nmin)
{
	facen = std::vector<std::vector<Face>>(4);
}



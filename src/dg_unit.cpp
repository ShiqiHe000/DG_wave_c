#include "dg_unit.h"
#include "dg_param.h"
#include <vector>


/// @brief
/// constructor (default).
/// default: level = 0, n = min poly order, m = max poly order
/// index = {0, 0, 0}, mpi_f = false, x/ycoords = {0.0, 0.0}
/// solution = nullptr. 
Unit::Unit() : n(grid::nmin), m(grid::nmin)
{
	facen = std::vector<std::vector<Face>>(4);

	ref_x = std::vector<double>(2);
	ref_y = std::vector<double>(2);

	// initialize element reference boudaries. Reference space [-1, 1]
	ref_x[0] = -1.0; ref_y[0] = -1.0;
	ref_x[1] =  1.0; ref_y[1] =  1.0;

}



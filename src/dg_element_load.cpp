#include "dg_param.h"
#include "dg_element_load.h"
#include <cmath>
#include <cassert>

/// @brief 
/// Element computational load. The load on each element due to fluid computations is O(N**4),
/// where N is the number of grid points along one direction. The load is normalized between 1 and 2. 
/// @param porder input the polynomial order of current element. (Now we assume polynomial
/// order are identical on two direction.)
double Elem_load(int porder){

	static const double load_min = std::pow((double)(grid::nmin + 1), 4);
	static const double load_max = std::pow((double)(grid::nmax + 1), 4);

	assert(load_min <= load_max && "Error: min poly order should be no less than max poly order. \n");
		
	double load{};

	if(load_min < load_max){
		load = (std::pow(porder + 1, 4) - load_min) / (load_max - load_min) + 1.0;
	}
	else{	// load_min == load_max, i.e., max poly order == min poly order

		load = 1.0;
	}


	return load;

}

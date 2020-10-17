#include "dg_CFL.h"
#include <cassert>
#include "dg_param.h"
#include <iostream>

/// @brief
/// Check CFL condition. t <= 1.73 * C / N ^ 2.  
/// C = 10 for Legendre polynomial. 
/// Polynomial order
void Check_CFL(){


	double del_t = 1.73 * 10 / (grid::nmax * grid::nmax);

	double t_user = dg_time::t_total / dg_time::nt;

	if(del_t < t_user){

		std::cout << "time step should be less than " << del_t << "\n";;
	}

	assert(del_t >= t_user && "Time step too large. Abort. \n");
}

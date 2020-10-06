#include "dg_CFL.h"
#include <cassert>
#include "dg_param.h"
#include <cmath>
#include <iostream>

void Check_CFL(){

	double d_length = grid::gx_r - grid::gx_l;

	int exp = grid::exp_x + grid::hlevel_max;

	double del_t = 1.73 * d_length / (std::pow(2, exp) * std::pow(grid::nmax, 2) * dg_fun::C);

	double t_user = dg_time::t_total / dg_time::nt;

	if(del_t < t_user){

		std::cout << "time step should be less than " << del_t << "\n";;
	}

	assert(del_t >= t_user && "Time step too large. Abort. \n");
}

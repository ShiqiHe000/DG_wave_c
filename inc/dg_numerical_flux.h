#ifndef DG_NUMERICAL_FLUX_H
#define DG_NUMERICAL_FLUX_H

#include <vector>

void Riemann_solver_x(std::vector<double>& q_l, std::vector<double>& q_r, 
			std::vector<double>& n_flux, int normal, std::vector<int>& index);
#endif

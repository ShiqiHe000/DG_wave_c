#ifndef DG_RIEMANN_SOLVER_H
#define DG_RIEMENN_SOLVER_H

#include <vector>

void Riemann_solver_x(std::vector<double>& q_l, std::vector<double>& q_r, 
			std::vector<double>& n_flux, int normal, std::vector<int>& index);

void Riemann_solver_y(std::vector<double>& q_l, std::vector<double>& q_r, 
			std::vector<double>& n_flux, int normal, std::vector<int>& index);
#endif

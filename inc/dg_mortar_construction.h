#ifndef DG_MORTAR_CONSTRUCTION_H
#define DG_MORTAR_CONSTRUCTION_H

#include <vector>

void L2_projection_to_mortar(int J, int n, int level, int l_max, double a, double b,
			 	std::vector<double>& solution_int, std::vector<double>& psi, 
				std::vector<double>& mapped_points);

void L2_projection_to_element(int J, int n, int level, int l_max, double a, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar, 
				std::vector<double>& mapped_points);
#endif

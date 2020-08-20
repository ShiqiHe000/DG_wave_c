#ifndef DG_SPATIAL_DERIVATIVE_H
#define DG_SPATIAL_DERIVATIVE_H

#include <vector>
#include <unordered_map>
#include "dg_unit.h"

void Spatial_derivative(int porder, std::unordered_map<int, std::vector<double>>& flux, 
			std::unordered_map<int, std::vector<double>>& flux_der, 
			Unit* temp, std::vector<int>& index);
#endif

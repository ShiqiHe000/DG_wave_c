#ifndef DG_FLUX_VECTOR_H
#define DG_FLUX_VECTOR_H

#include <vector>
#include <unordered_map>

void xflux(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& xf, int index);

void yflux(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& yf, int index);

#endif

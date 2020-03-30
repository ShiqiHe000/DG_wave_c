#include "dg_flux_vector.h"
#include <vector>
#include "dg_param.h"

/// @brief
/// Compute the horizontal fluxes. X direction. 
/// @param q solutions.
/// @param xf horizontal numerical fluxes.
/// @parma solution indeices (3 equaitons).
void xflux(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& xf, int index){

	xf[0].push_back(dg_fun::C * dg_fun::C * q[1][index]);

	xf[1].push_back(q[0][index]);

	xf[2].push_back(0.0);

}

/// @brief
/// Compute the horizontal fluxes. Y direction. 
/// @param q solutions.
/// @param xf horizontal numerical fluxes.
/// @parma solution indeices (3 equaitons).
void yflux(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& yf, int index){

	yf[0].push_back(dg_fun::C * dg_fun::C * q[2][index]);

	yf[1].push_back(0.0);

	yf[2].push_back(q[0][index]);

}

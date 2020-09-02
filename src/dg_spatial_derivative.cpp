#include "dg_spatial_derivative.h"
#include "dg_unit.h"
#include "dg_param.h"
#include <vector>
#include "dg_basis.h"
#include "dg_nodal_2d_storage.h"
#include <unordered_map>
#include <iostream>	// test
#include "dg_single_index.h"	// test


/// @brief
/// Compute the first spatial derivative given by the interior points and two boundary points.
/// @param flux Flux.
/// @param flux_der flux derivatives
/// @param temp pointer to the current element.
/// @param index node index for boundary numerical fluxes. 
void Spatial_derivative(int porder, std::unordered_map<int, std::vector<double>>& flux, 
			std::unordered_map<int, std::vector<double>>& flux_der, 
			Unit* temp, std::vector<int>& index){
	// compute flux derivatives
	for(int s = 0; s < dg_fun::num_of_equation; ++s){

		flux_der[s] = std::vector<double>(porder + 1);
		Matrix_vector_multiplication(porder, nodal::mfirst_der[porder], flux[s], flux_der[s]);

	}


	for(int s = 0; s < dg_fun::num_of_equation; ++s){

		for(int j = 0; j <= porder; ++j){

			flux_der[s][j] += ((temp -> nflux_r[index[s]]) * nodal::lagrange_r[porder][j] 
						+ (temp -> nflux_l[index[s]]) * nodal::lagrange_l[porder][j]) 
						/ nodal::gl_weights[porder][j];

		}
	
	}

	
	
}

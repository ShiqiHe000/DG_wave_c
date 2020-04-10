#include "dg_unit.h"
#include "dg_mortar_construction.h"
#include <vector>
#include "dg_param.h"
#include "dg_basis.h"
#include <cmath>	// pow
#include <algorithm>
#include "dg_nodal_2d_storage.h"

/// @brief
/// L2 projection from element to mortar.
/// @param J Maximum polynomial order between left and right element. 
/// @parma n Element polynomial order.
/// @param level element h-refinement level.
/// @param l_max Maximum h-refinement level between two elements. 
/// @parma a Mortar coorindate offset. 
/// @parma b Mortar coorindate scaling. 
/// @parma solution_int Element interface solution. 
/// @param psi Mortar solution. 
void L2_projection_to_mortar(int J, int n, int level, int l_max, double a, double b,
			 	std::vector<double>& solution_int, std::vector<double>& psi){

	if(J == n && level == l_max){	// direct copy

		// U = psi
		psi = solution_int;	
	}
	else{

		std::vector<int> index_elem{0, n + 1, (n + 1) * 2}; // 3 equation
		std::vector<int> index_mortar{0, J + 1, (J + 1) * 2}; // 3 equation

		std::vector<double> bary(n + 1);

		BARW(n, nodal::gl_points[n], bary);
		
		// L2 projection
		for(int i = 0; i <= J; ++i){	

			double s = a + b * nodal::gl_points[J][i];	// map GL point from mortar to element

			std::vector<double> lag(n + 1);

			// get the lagrange interpolation value at this point (s).
			Lagrange_interpolating_polynomial(n, s, nodal::gl_points[n], bary, lag);

			for(int j = 0; j <= n; ++j){

				double inter = lag[j] * std::pow(nodal::gl_weights[J][i], 2);

				for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

					psi[index_mortar[equ]] += inter * solution_int[index_elem[equ]];
				}

				std::transform(index_elem.begin(), index_elem.end(), index_elem.begin(), 
						[](int x){return x + 1;});
			}

			std::transform(index_mortar.begin(), index_mortar.end(), index_mortar.begin(), 
					[](int x){return x + 1;});
		}


	}


}



/// @brief
/// L2 projection from element to mortar.
/// @param J Maximum polynomial order between left and right element. 
/// @parma n Element polynomial order.
/// @param level element h-refinement level.
/// @param l_max Maximum h-refinement level between two elements. 
/// @parma a Mortar coorindate offset. 
/// @parma b Mortar coorindate scaling. 
/// @parma solution_int Element interface solution. 
/// @param psi Mortar solution. 
void L2_projection_to_element(int J, int n, int level, int l_max, double a, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar){

	if(J == n && level == l_max){	// direct copy

		// U = psi
		nflux_elem = nflux_mortar;	
	}
	else{

		std::vector<int> index_elem{0, n + 1, (n + 1) * 2}; // 3 equation
		std::vector<int> index_mortar{0, J + 1, (J + 1) * 2}; // 3 equation

		std::vector<double> bary(J + 1);

		BARW(n, nodal::gl_points[J], bary);
		
		// L2 projection
		for(int i = 0; i <= n; ++i){	

			double z = (nodal::gl_points[n][i] - a) / b;	// map GL point from elem to mortar

			std::vector<double> lag(J + 1);

			// get the lagrange interpolation value at this point (s).
			Lagrange_interpolating_polynomial(J, z, nodal::gl_points[J], bary, lag);

			for(int j = 0; j <= J; ++j){

				double inter = lag[j] * std::pow(nodal::gl_weights[n][i], 2);

				for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

					nflux_elem[index_elem[equ]] += inter * nflux_mortar[index_mortar[equ]];
				}

				std::transform(index_mortar.begin(), index_mortar.end(), index_mortar.begin(), 
						[](int x){return x + 1;});
			}

			std::transform(index_elem.begin(), index_elem.end(), index_elem.begin(), 
					[](int x){return x + 1;});
		}


	}


}

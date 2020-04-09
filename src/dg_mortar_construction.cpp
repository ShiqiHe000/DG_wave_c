#include "dg_unit.h"
#include "dg_mortar_construction.h"
#include <vector>
#include "dg_param.h"
#include "dg_basis.h"
#include <cmath>	// pow
#include <algorithm>
#include "dg_nodal_2d_storage.h"


/// @param temp pointer to the current element
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
		
		// generate the projection matrix
		for(int i = 0; i <= J; ++i){	

			double inter{};

			double s = a + b * nodal::gl_points[J][i];	// map GL point from mortar rto element

			std::vector<double> lag(n + 1);

			// get the lagrange interpolation value at this point (s).
			Lagrange_interpolating_polynomial(n, s, nodal::gl_points[n], bary, lag);

			for(int j = 0; j <= n; ++j){

				inter = lag[j] * std::pow(nodal::gl_weights[J][i], 2);

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

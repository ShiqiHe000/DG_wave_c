#include "dg_error_estimator.h"
#include "dg_unit.h"
#include <vector>
#include <cassert>
#include "dg_param.h"
#include "dg_single_index.h"
#include <unordered_map>
#include "dg_basis.h"
#include "dg_nodal_2d_storage.h"

/// @brief
/// 
/// @param temp pointer to the estimated element.
/// @note The polynomial order in x and y direction should be identical.
double Error_indicator(Unit* temp){

	assert((temp -> n > 6) && "Polynomial order is too low to estimate error.");

	std::unordered_map<int, std::vector<double>> ap;	// discrete spectrum ap
	for(int i = 0; i < dg_fun::num_of_equation; ++i){	// allocate space

		ap[i] = std::vector<double> (dg_refine::fit_point_num);
	}


	int N = temp -> n;
	int M = temp -> m;
	int p = N;	// assume N == M

	for(int node = dg_refine::fit_point_num - 1; node >= 0; --node){

		// x direction-----------------------------------------------------
		for(int i = 0; i <= p; ++i ){	// compute sum of coefficients: sum(a(i, p))
			
			// a(n, m) = a(i, p), apply Gauss quadrature
			for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

				for(int k = 0; k <= N; ++k){	// x dir
	
					// x
					// compute Legendre polynomial of order i at GL point k
					double legendre_x{};
					double legendre_dir{};
					Legendre_polynomial_and_derivative(i, nodal::gl_points[N][k], 
										legendre_x, legendre_dir);
					for(int l = 0; l <= M; ++l){	// y dir
	
						int index = Get_single_index(k, l, M + 1);

						double legendre_y{};
						
						// y
						// compute Legendre polynomial of order l at GL point l
						Legendre_polynomial_and_derivative(p, nodal::gl_points[M][l], 
										legendre_y, legendre_dir);

						ap[equ][node] += (temp -> solution[equ][index]) * 
								legendre_x * legendre_y * 
								nodal::gl_weights[N][k] * nodal::gl_weights[M][l];
								
						
					}
					
				}

				ap[equ][node] *= (2.0 * (double)i + 1.0) * (2.0 * (double)p + 1.0) / 4.0;
			}

		}
		//-----------------------------------------------------------------

		// y direction-----------------------------------------------------
		for(int j = 0; j < p; ++j ){	// here does not need to include the last point (j < p)
			
					
			for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

				for(int k = 0; k <= N; ++k){	// x dir
	
					// x
					// compute Legendre polynomial of order i at GL point k
					double legendre_x{};
					double legendre_dir{};
					Legendre_polynomial_and_derivative(p, nodal::gl_points[N][k], 
										legendre_x, legendre_dir);
					for(int l = 0; l <= M; ++l){	// y dir
	
						int index = Get_single_index(k, l, M + 1);

						double legendre_y{};
						
						// y
						// compute Legendre polynomial of order l at GL point l
						Legendre_polynomial_and_derivative(j, nodal::gl_points[M][l], 
										legendre_y, legendre_dir);

						ap[equ][node] += (temp -> solution[equ][index]) * 
								legendre_x * legendre_y * 
								nodal::gl_weights[N][k] * nodal::gl_weights[M][l];
								
						
					}
					
				}

				ap[equ][node] *= (2.0 * (double)j + 1.0) * (2.0 * (double)p + 1.0) / 4.0;

			}

		}
		//-----------------------------------------------------------------
		p -= 2;
	}
		
}

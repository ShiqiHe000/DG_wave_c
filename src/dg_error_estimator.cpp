#include "dg_error_estimator.h"
#include "dg_unit.h"
#include <vector>
#include <cassert>
#include "dg_param.h"
#include "dg_single_index.h"
#include <unordered_map>
#include "dg_basis.h"
#include "dg_nodal_2d_storage.h"
#include <cmath>
#include "dg_local_storage.h"
#include <iostream>	// test

// forward declaration---------------------------------------------------------
double Decay_rate(std::vector<int>& porder, std::vector<double>& ap, double& C);

void Error_indicator(Unit* temp, std::vector<double>& sigma, std::vector<bool>& flag_refine, std::vector<bool>& flag_coarsen);

double Solution_l2_norm(int equ, Unit* temp);

void Refinement_flag();
//-----------------------------------------------------------------------------

void Refinement_flag(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		std::vector<double> sigma(dg_fun::num_of_equation);
		std::vector<bool> flag_refine(dg_fun::num_of_equation);
		std::vector<bool> flag_coarsen(dg_fun::num_of_equation);

		Error_indicator(temp, sigma, flag_refine, flag_coarsen);

		// now refine depends on pressure--------------
		if(flag_refine.front()){	// refine depends on pressure

			if(sigma.front() < 1){	// h-refinemnt
			
				if(temp -> index[2] < grid::hlevel_max){ // if no exceed the max hlevel
					temp -> hrefine = true;
				}
	
			}
			else{	// p-refinement
				if(temp -> n < grid::nmax){	// not exceed the max poly order
					temp -> prefine = true;
				}
			}
		}
		else if(flag_coarsen.front()){

			temp -> coarsen = true;	// coarsening
			// test no coarsening 
//			temp -> coarsen = false;	

		}
		//-------------------------------------------------
		temp = temp -> next;
	}

}

/// @brief
/// Calculate the error indicator sigma, if sigma < 1 apply h-refinemnt, 
/// if sigma >= 1 apply prefinement. 
/// @param temp pointer to the estimated element.
/// @param sigma the error indicator of each function. 
/// @param flag_refine refinement flag.
/// @param flag_coarsen coarsening flag.
/// @note The polynomial order in x and y direction should be identical.
void Error_indicator(Unit* temp, std::vector<double>& sigma, std::vector<bool>& flag_refine, std::vector<bool>& flag_coarsen){

	assert((temp -> n >= 4) && "Polynomial order is too low to estimate error.");

	std::unordered_map<int, std::vector<double>> ap;	// discrete spectrum ap
	for(int i = 0; i < dg_fun::num_of_equation; ++i){	// allocate space

		ap[i] = std::vector<double> (dg_refine::fit_point_num);
	}


	int N = temp -> n;
	int M = temp -> m;
	int p = N;	// assume N == M

	assert((M == N) && "Right now the error estimator only works for same order in 2D. ");

	std::vector<int> porder(dg_refine::fit_point_num);	// record the poilynomial order

	for(int node = dg_refine::fit_point_num - 1; node >= 0; --node){

		porder[node] = p;

		// x direction-----------------------------------------------------
		for(int i = 0; i <= p; ++i ){	// compute sum of coefficients: sum(a(i, p))
			
			// a(n, m) = a(i, p), apply Gauss quadrature
			for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

				double ap_now{};

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

						ap_now += (2.0 * (double)i + 1.0) * (2.0 * (double)p + 1.0) / 4.0 *
								(temp -> solution[equ][index]) * 
								legendre_x * legendre_y * 
								nodal::gl_weights[N][k] * nodal::gl_weights[M][l];
					}
					
				}
	
				ap[equ][node] += std::abs(ap_now);

			}

		}
		//-----------------------------------------------------------------

		// y direction-----------------------------------------------------
		for(int j = 0; j < p; ++j ){	// here does not need to include the last point (j < p)
			
					
			for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

				double ap_now{};

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

						ap_now += (2.0 * (double)j + 1.0) * (2.0 * (double)p + 1.0) / 4.0 *
								(temp -> solution[equ][index]) * 
								legendre_x * legendre_y * 
								nodal::gl_weights[N][k] * nodal::gl_weights[M][l];
						
					}
					
				}

				ap[equ][node] += std::abs(ap_now);
			}

		}
		//-----------------------------------------------------------------
		--p;
	}
		

	// refinment criteria: if exceed the acceptable level then flag as needed refinement. 
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		// L2 norm of the current element =====================================================
		double u_norm = Solution_l2_norm(equ, temp);

		double tol_min = u_norm * dg_refine::tolerance_min;	// the threshold for refinement 
		
		double tol_max = u_norm * dg_refine::tolerance_max;	// the threshold for coarening		
		// ====================================================================================


		// use threshold as the refinement criterion =========================================
//
//		double tol_min = dg_refine::tolerance_min;	// the threshold for refinement 
//		
//		double tol_max = dg_refine::tolerance_max;	// the threshold for coarening		

		// ===================================================================================

		double C{};	// the const: C*exp(-sigma * n)

		sigma[equ] = Decay_rate(porder, ap[equ], C);
		
		// sum of error
		double sum = std::sqrt(std::pow(ap[equ].back(), 2) 
				+ C * C / (2.0 * sigma[equ]) * std::exp(- 2.0 * sigma[equ] * (porder.back() + 1)));
		
//		double sum = ap[equ].back();	// sum of the last spectrum

		if(sum > tol_min){	// need refine

			flag_refine[equ] = true;

			flag_coarsen[equ] = false;

//			// get decay indicator
//			sigma[equ] = Decay_rate(porder, ap[equ]);

		}
		else if(sum <= tol_max ){	// need coarsen

			flag_refine[equ] = false;

			flag_coarsen[equ] = true;
			

		}
		else{	// if error in between then do nothing

			flag_refine[equ] = false;

			flag_coarsen[equ] = false;
		}


	}
	

}

/// @brief
/// Calculate the L2 norm of solution of the current element of the demanding equaiton.
/// @param equ the equation number
/// @param temp pointer to the Current element.
double Solution_l2_norm(int equ, Unit* temp){

	double norm{};

	int nodei{};

	for(int i = 0; i <= temp -> n; ++i){

		double weight_x = nodal::gl_weights[temp -> n][i];

		for(int j = 0; j <= temp -> m; ++j){

			norm += std::pow(temp -> solution[equ][nodei], 2) * nodal::gl_weights[temp -> m][j] * weight_x;
			++nodei;
		}

	}


	return (std::sqrt(norm));

}



/// @brief
/// Assume the decaying spectrum ap is an exponential function ( ap  = const * exp ^ (-sigma * n)).
/// Apply log arithmetic to the two side of the function and fit the ap by a linear regression line.
/// Then we get the decay factor (the absolute value of the line). 
/// @param porder the corresponding order of ap. (x)
/// @param ap the decaying spectrums. (y)
/// @param C the constant in the exponential function. 
double Decay_rate(std::vector<int>& porder, std::vector<double>& ap, double& C){

	double x_avg{}; double y_avg{};

	for(int i = 0; i < dg_refine::fit_point_num; ++i){

		x_avg += (double)porder[i];
		y_avg += std::log(ap[i]);
	}
	x_avg /= (double)dg_refine::fit_point_num;
	y_avg /= (double)dg_refine::fit_point_num;

	double sigma{};
	double numer{};
	double denumer{};

	for(int i = 0; i < dg_refine::fit_point_num; ++i){
		
		numer += ((double)porder[i] - x_avg) * (std::log(ap[i]) - y_avg);
		denumer += std::pow(((double)porder[i] - x_avg), 2);
	
	}
		
	sigma = numer / denumer;

	double alpha = y_avg - sigma * x_avg; // alpha = ln(C)

	C = std::exp(alpha);

	return std::abs(sigma);

}

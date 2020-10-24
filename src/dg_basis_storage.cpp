#include <vector>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_basis.h"
#include "dg_basis_storage.h"
#include "dg_single_index.h"
#include "dg_write_first_der_matrix.h"
#include "dg_eigenvalues.h"
#include <iostream>	// test


/// @brief
/// Processor's local storages.
/// Generate the basis information (GL PIOINTS, WEIGHTS, DERIVATIVE 
/// MATRIX, etc.) for all polynomial level at the beginning,
/// and store it in hash tables. 
/// Therefore, after refinement start we do not need to 
/// re-generate the basis information, or store duplicated data. 
void Construct_basis_storage(){


	// generate the 2d storages
	for(int n = grid::nmin; n <= grid::nmax; n+=2){	// k = poly order
		
		// allocate space	
		nodal::gl_points[n] = std::vector<double>(n + 1);
		nodal::gl_weights[n] = std::vector<double>(n + 1);
	
		int der_size = (n + 1) * (n + 1);
		nodal::first_der[n] = std::vector<double>(der_size);
		nodal::mfirst_der[n] = std::vector<double>(der_size);
	
		// generate current gl_p, gl_w
		GL(n, nodal::gl_points[n], nodal::gl_weights[n]);
	
		// compute GLL points, weights--------------------------------
//		nodal::gll_points[n] = std::vector<double>(n + 1);
//		nodal::gll_weights[n] = std::vector<double>(n + 1);
//		GLL(n, nodal::gll_points[n], nodal::gll_weights[n]);
		//------------------------------------------------------------
	
		// compute Chebyshev GLL points, weights--------------------------------
		nodal::cgll_points[n] = std::vector<double>(n + 1);
		nodal::cgll_weights[n] = std::vector<double>(n + 1);
		nodal::chebyshev_first_der[n] = std::vector<double>(der_size);

//		Chebyshev_GLL_nodes_and_weights(n, nodal::cgll_points[n], nodal::cgll_weights[n]);
//		Chebyshev_GLL_first_der_matrix(n, nodal::cgll_points[n], nodal::chebyshev_first_der[n]);
	
		Chebyshev_GLL_nodes_Baltensperger(n, nodal::cgll_points[n]);
		Chebyshev_first_der_trigonometric(n, nodal::cgll_points[n], nodal::chebyshev_first_der[n]);

		Eigen_Chebyshev_output(n);	// compute eigenvalues and output to files. 
		//------------------------------------------------------------

		std::vector<double> bary(n + 1);
	
		BARW(n, nodal::gl_points[n], bary);
	
		// first order derivative matrix ====================================================================
		//Mth_order_polynomial_derivative_matrix(n, 1, nodal::gl_points[n], nodal::first_der[n], bary);

		// minimize the round-off error version
		First_order_polynomial_derivative_matrix(n, nodal::gl_points[n], nodal::first_der[n], bary);
		// ==================================================================================================

		// compute eigenvalues of the first derivative matrix and output to files -------------------------
		//Compute_eigs_plus_output(n);
		//-------------------------------------------------------------------------------------------------

		// output the first der matrix to file --------------------
		//Write_first_der_matrix(n);
		// -------------------------------------------------------

		// Modify first derivative
		for(int j = 0; j <= n; ++j){

			for(int i = 0; i <= n; ++i){

				int index1 = Get_single_index(i, j, n + 1);
				int index2 = Get_single_index(j, i, n + 1);

				nodal::mfirst_der[n][index1] = - nodal::first_der[n][index2] * nodal::gl_weights[n][j] / 
								nodal::gl_weights[n][i];

			}

		}
		// ---------------------------------
		//Write_mfirst_der_matrix(n);
		// ---------------------------------
		
		// DG first derivative matrix---------
	//	Compute_DG_first_der_eigs_plus_output(n);
		// ------------------------------------


		// first der matrix is not needed anymore.
		nodal::first_der.clear();

		// Lagrange interpolates on the boundaries. 
		nodal::lagrange_l[n] = std::vector<double>(n + 1);
		nodal::lagrange_r[n] = std::vector<double>(n + 1);
		Lagrange_interpolating_polynomial(n, -1.0, nodal::gl_points[n], bary, nodal::lagrange_l[n]);
		Lagrange_interpolating_polynomial(n,  1.0, nodal::gl_points[n], bary, nodal::lagrange_r[n]);

	}

	// if we apply refinement, we need the value of Legendre polynomial 
//	if(dg_refine::adapt){
//
//		for(int n = grid::nmin; n <= grid::nmax; n+=2){	// k = poly order
//			
//			nodal::legendre[n] = std::vector<double>(n + 1);	// legendre polynomial for error estimator
//
//			// Legendre polynomial value at GL points
//			for(int i = 0; i <= n; ++i ){
//				double dq{}; 	// legendre poly derivative, useless here.
//				Legendre_polynomial_and_derivative(n, nodal::gl_points[n][i], nodal::legendre[n][i], dq);
//			}
//		}
//
//
//		int n = grid::nmin - 2;
//
//		for(; n >= 1; n-=2 ){
//
//			// allocate space	
//			nodal::gl_points[n] = std::vector<double>(n + 1);
//			nodal::gl_weights[n] = std::vector<double>(n + 1);
//
//			nodal::legendre[n] = std::vector<double>(n + 1);	// legendre polynomial for error estimator
//
//			// generate current gl_p, gl_w
//			GL(n, nodal::gl_points[n], nodal::gl_weights[n]);
//
//			for(int i = 0; i <= n; ++i ){
//				double dq{}; 	// legendre poly derivative, useless here.
//				Legendre_polynomial_and_derivative(n, nodal::gl_points[n][i], nodal::legendre[n][i], dq);
//			}
//		}
//
//	}


}


#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "dg_single_index.h"
#include "dg_eigenvalues.h"

/// @brief
/// Compute the eigenvalues of the first derivative matrix. 
/// @param n polynomial order.
/// @param der first derivative matrix.
/// @param eigenvalues results. 
void First_der_eigens(int n, std::vector<double>& der, 
			std::vector<double>& eig_real, std::vector<double>& eig_im){

	Eigen::MatrixXd m(n + 1, n + 1);

	// define first der matrix for Eigen. 
	for(int i = 0; i <= n; ++i){

		for(int j = 0; j <= n; ++j){

			int nodei = Get_single_index(i, j, n + 1);

			m(i, j) = der[nodei];
		}
	}

	// Solve
	Eigen::EigenSolver<Eigen::MatrixXd> eig_res(m);

	// Check valid
	if(eig_res.info() != Eigen::Success) std::abort();

	// transfer into vector
	for(int i = 0; i <= n; ++i){

		eig_real.push_back(eig_res.eigenvalues()[i].real());	// get the real part.
		eig_im.push_back(eig_res.eigenvalues()[i].imag());	// get the imaginary part. 
	}

}

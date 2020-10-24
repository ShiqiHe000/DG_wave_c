#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "dg_single_index.h"
#include "dg_eigenvalues.h"
#include <mpi.h>
#include "dg_param.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "dg_nodal_2d_storage.h"

// formward declaration --------------------------------------------------------------
void Compute_eigs_plus_output(int n);

void Compute_DG_first_der_eigs_plus_output(int n);

void Eigen_Chebyshev_output(int n);

void First_der_eigens(int n, std::vector<double>& der, 
			std::vector<double>& eig_real, std::vector<double>& eig_im);

void Write_eigs(int n, std::vector<double>& eig_real, std::vector<double>& eig_im);
// -----------------------------------------------------------------------------------

/// @brief
/// Compute the standard first derivative's eigenvalues and output to files. 
/// @param n polynomial order.
void Compute_eigs_plus_output(int n){

	std::vector<double> eig_real;
	std::vector<double> eig_im;

	First_der_eigens(n, nodal::first_der[n], eig_real, eig_im);

	Write_eigs(n, eig_real, eig_im);
}

/// @brief
/// Compute the DG first derivative's eigenvalues and output to files. 
/// @param n polynomial order.
void Compute_DG_first_der_eigs_plus_output(int n){

	std::vector<double> eig_real;
	std::vector<double> eig_im;

	First_der_eigens(n, nodal::mfirst_der[n], eig_real, eig_im);

	Write_eigs(n, eig_real, eig_im);
}

/// @brief
/// Compute the Chebyshev first derivative's eigenvalues and output to files. 
/// @param n polynomial order.
void Eigen_Chebyshev_output(int n){

	std::vector<double> eig_real;
	std::vector<double> eig_im;

	First_der_eigens(n, nodal::chebyshev_first_der[n], eig_real, eig_im);

	Write_eigs(n, eig_real, eig_im);
}


/// @brief
/// Compute the eigenvalues of the first derivative matrix. 
/// @param n polynomial order.
/// @param der first derivative matrix.C
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

	// matrix transpose
//	m.transposeInPlace();

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

/// @brief
/// Write eigenvalues to file
/// @param n polynomial order.
/// @param eig_real real parts of the eigenvalues. 
/// @param eig_im imaginary parts of the eigenvalues. 
void Write_eigs(int n, std::vector<double>& eig_real, std::vector<double>& eig_im){


	if(mpi::rank == 0){

		// generate the file name
		std::stringstream ss;
		ss << fileinfo::eigenvalues_place <<"eig" << std::setfill('0') << std::setw(2) << n << ".txt";
		std::string filename = 	ss.str();
		std::ofstream myfile; 	// stream class to write on files	

		// processor open file
		myfile.open(filename, std::ios::trunc);	// truncate the old file

		myfile<< std::fixed;
		myfile<< std::setprecision(17);

		int index{};

		for(int i = 0; i <= n; ++i){

			myfile << eig_real[i] << " " << eig_im[i] << "\n";

		}

		myfile.close();
		

	}
}

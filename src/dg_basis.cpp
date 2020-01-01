#include <iostream>
#include <limits>	// epsilon
#include <cmath>	// sqrt
#include "dg_single_index.h"


const double pi = 4.0 * atan(1.0); 

/// @brief 
/// compute matrix vector multiplication.
/// d * f = der. Algorithm 19.
/// @param n polynomial order
/// @param d coefficient martrix (usually derivative matrix), d[n * n].
/// @param f vector, f[n].
/// @param der output. (usually the derivative of interpolation), der[n].
void Matrix_vector_multiplication(int n, double* d, double* f, double* der ){

//	for(){
//	
//
//
//	}

}



/// @brief 
/// compute matrix vector multiplication.
/// d * f = der. Use BLAS library. 
/// @param n polynomial order
/// @param d coefficient martrix (usually derivative matrix), d[n * n].
/// @param f vector, f[n].
/// @param der output. (usually the derivative of interpolation), der[n].
void Matrix_vector_multiplication_blas(int n, double* d, double* f, double* der){


}


/// @brief
/// Legendre polynomial of degree k and its derivative using the three
/// term recursive. 
/// Algorithm 22
/// @param n polynomial order
/// @param x GL points
/// @param q Legendre polynomial of degree k
/// @param dq Derivative of Legendre polynomial  
void Legendre_polynomial_and_derivative(int n, double& x, double& q, double& dq){
	if(n == 0){
		q = 1.0;
		dq = 0.0;
	}
	else if(n == 1){
		q = x;
		dq = 1.0;

	}
	else{
		double q_m2 = 1.0;
		double q_m1 = x;
		double dq_m2 = 0.0;
		double dq_m1 = 1.0;
		
		for(int k = 2; k <= n; ++k ){
			q = (double)(2 * k - 1) * x * q_m1 / (double)(k) - (double)(k -1) * q_m2 / (double)(k);
			dq = dq_m2 + (double)(2 * k - 1) * q_m1;
			q_m2 = q_m1;
			q_m1 = q;
			dq_m2 = dq_m1;
			dq_m1 = dq;
		}

	}
	
}


/// @brief
/// Compute the Gauss Legendre nodes and weights. Algprithm 23
/// @param n polynomial order
/// @param gl_p GL points
/// @param gl_w GL weigths
void GL(int n, double* gl_p, double* gl_w){

	double delta;
	double q, dq, tol;

	tol = 4.0 * std::numeric_limits<double>::epsilon();
	q = 0.0;
	dq = 0.0;
	
	if(n == 0){
		gl_p[0] = 0.0;
		gl_w[0] = 2.0;
	}
	else if(n == 1){

		gl_p[0] = - sqrt(1.0 / 3.0);
		gl_w[0] = 1.0;
		gl_p[1] = sqrt(1.0 / 3.0);
		gl_w[1] = 1.0;

	}
	else{
		for(int j = 0; j <= ((n+1)/2)-1; ++j ){
			// initial guess
			gl_p[j] = - cos(pi * (double)(2 * j + 1)/(double)(2 * n + 2));

			// iterative method
			delta = 10000000.0;
			while(true){
				
				Legendre_polynomial_and_derivative(n+1, gl_p[j], q, dq);

				delta = - q / dq;
				gl_p[j] = gl_p[j] + delta;
				
				if(std::abs(delta) <= (tol*std::abs(gl_p[j])))
					break;
			}

			Legendre_polynomial_and_derivative(n+1, gl_p[j], q, dq);
				
			gl_p[n - j] = - gl_p[j];
			gl_w[j] = 2.0 / ((1.0 - pow(gl_p[j], 2.0)) * pow(dq, 2.0));
			gl_w[n - j] = gl_w[j];

		}

	}


	if(n % 2 == 0){
		double gl_p_now = 0.0;
		Legendre_polynomial_and_derivative(n+1, gl_p_now, q, dq);

		gl_p[n/2] = 0.0;
		gl_w[n/2] = 2.0 / pow(dq, 2.0);

	}


}

/// @brief 
/// Barycentric weights for Lagrange Iterpolation. Algorithm 30.
/// @param n polynomial order
/// @param x spectral points
/// @param bary barycentric weights
void BARW(int n, double* x, double* bary){
	
	for(int i = 0; i <= n; ++i){
		bary[i] = 1.0;
	}

	for(int j = 1; j <= n; ++j){
		for(int k = 0; k <= j-1; ++k){
			bary[k] = bary[k] * (x[k] - x[j]);
			bary[j] = bary[j] * (x[j] - x[k]);

		}

	}

	for(int j = 0; j <= n; ++j){

		bary[j] = 1.0 / bary[j];

	}


}



/// @brief
/// M-th order derivative matrix. Algorithm 36-37
/// @param n polynomial order
/// @param mth_der m-th order polynomial derivative
/// @param x spectral points
/// @param der m-th order derivative matrix
void Mth_order_polynomial_derivative_matrix(int n, int mth_der, double* x, double* der){
	
	double* bary = new double[n+1];

	BARW(n, x, bary);

	// mth-order == 1
	for(int i = 0; i <= n; ++i){
		int inode = Get_single_index(i, i, n+1);
		der[inode] = 0.0;
		
		for(int j = 0; j <= n; ++j){
			if(j != i){
				int node1 = Get_single_index(i, j, n+1);
				int node2 = Get_single_index(i, i, n+1);

				der[node1] = bary[j] / bary[i] / (x[i] - x[j]);
				der[node2] = der[node2] - der[node1];

			}

		}

	}	
	

	if(mth_der == 1){
		return;
	}


	// mth_order > 1
	//------------
//	aux = der
//	//------------
//	for(int k = 2; k <= mth_der; ++k){
//		for(int i = 0; i <= n;, ++i){
//			
//			int inode = Get_single_index(i, i, n+1);
//			der[inode] = 0.0;
//
//			for(int j = 0; j <= n; ++j){
//				if(j != i){
//					
//					int node1 = Get_single_index(i, j, n+1);
//					int node2 = Get_single_index(i, i, n+1);
//
//					der[node1] = (double)(k) / (x[i] - x[j]) * 
//							(bary[j] / bary[i] * aux[node2] - aux[node1]);
//
//					der[node2] = der[node2] - der[node1];
//				}
//
//			}
//
//
//		}
//
//	}
	



}

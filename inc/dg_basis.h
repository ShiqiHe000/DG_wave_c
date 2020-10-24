#ifndef DG_BASIS_H
#define DG_BASIS_H

#include <vector>

bool Almost_equal(double a, double b);

void Mth_order_polynomial_derivative_matrix(int n, int mth_der, std::vector<double>& x, std::vector<double>& der, 
						std::vector<double>& bary);

void First_order_polynomial_derivative_matrix(int n, std::vector<double>& x, std::vector<double>& der, 
						std::vector<double>& bary);

void GL(int n, std::vector<double>& gl_p, std::vector<double>& gl_w);

void BARW(int n, std::vector<double>& x, std::vector<double>& bary);

double Lagrange_interpolation(int n, double x, std::vector<double>& xi, 
				std::vector<double>& f, std::vector<double>& w, std::vector<int>& index);

void Lagrange_interpolating_polynomial(int n, double target_p, std::vector<double>& x, std::vector<double>& bary,
					 std::vector<double>& lag );

double Interpolate_to_boundary(int n, std::vector<double>& q, std::vector<double>& lag);

void Matrix_vector_multiplication(int n, std::vector<double>& d, std::vector<double>& f, std::vector<double>& out);

void Legendre_polynomial_and_derivative(int n, double x, double& q, double& dq);

void GLL(int n, std::vector<double>& gll_p, std::vector<double>& gll_w);

void Chebyshev_GLL_nodes_and_weights(int n, std::vector<double>& cgll_p, std::vector<double>& cgll_w);

void Chebyshev_GLL_first_der_matrix(int n, std::vector<double>& cgll_p, std::vector<double>& cder);

void Chebyshev_first_der_trigonometric(int n, std::vector<double>& cgll_p, std::vector<double>& cder);

void Chebyshev_GLL_nodes_Baltensperger(int n, std::vector<double>& cgll_p);

#endif

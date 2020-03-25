#ifndef DG_BASIS_H
#define DG_BASIS_H

#include <vector>

void Mth_order_polynomial_derivative_matrix(int n, int mth_der, std::vector<double>& x, std::vector<double>& der, 
						std::vector<double>& bary);

void GL(int n, std::vector<double>& gl_p, std::vector<double>& gl_w);

void BARW(int n, std::vector<double>& x, std::vector<double>& bary);

void Lagrange_interpolating_polynomial(int n, double target_p, std::vector<double>& x, std::vector<double>& bary,
					 std::vector<double>& lag );


#endif

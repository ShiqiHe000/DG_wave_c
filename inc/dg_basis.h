#ifndef DG_BASIS_H
#define DG_BASIS_H

void GL(int n, double* gl_p, double* gl_w);


void BARW(int n, double* x, double* bary);

void Lagrange_interpolating_polynomial(int n, double target_p, double* x, double* bary, double* lag);

void Mth_order_polynomial_derivative_matrix(int n, int mth_der, double* x, double* der);



#endif

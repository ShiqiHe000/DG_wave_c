#include <mpi.h>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_basis.h"
#include "dg_poly_level_and_order.h"
#include "dg_basis_storage.h"
#include <iostream>	//test
// forward declaration----------------------------------------------------------------------------
void Get_nodal_2d_storage_basis(int n, int k, double** gl_p, double** gl_weight, double** first_der);

void Get_nodal_2d_storage_extends(int n, int plevel, double** lag_l, double** lag_r, double** gl_p);
//------------------------------------------------------------------------------------------------

/// @brief
/// Processor's local storages.
/// Generate the basis information (GL PIOINTS, WEIGHTS, DERIVATIVE 
/// MATRIX, etc.) for all polynomial level at the beginning,
/// and store it in hash tables. 
/// Therefore, after refinement start we do not need to 
/// re-generate the basis information, or store duplicated data. 
void Construct_basis_storage(){

	// generate gl_p and gl_w
	int level_max = Poly_order_to_level(grid::nmin, grid::nmax);

	nodal::gl_p = new double*[level_max + 1]{};	// level starts with 0
	nodal::gl_w = new double*[level_max + 1]{};

	// first derivative matrix
	nodal::first_der = new double*[level_max + 1]{};

	// lagrange interpolating polynomial
	nodal::lagrange_l = new double*[level_max + 1]{};
	nodal::lagrange_r = new double*[level_max + 1]{};

	// generate the 2d storages
	for(int k = 0; k <= level_max; ++k ){
		
		int porder = Poly_level_to_order(grid::nmin, k);
		
		Get_nodal_2d_storage_basis(porder, k, nodal::gl_p, nodal::gl_w, nodal::first_der);

		Get_nodal_2d_storage_extends(porder, k, nodal::lagrange_l, nodal::lagrange_r, nodal::gl_p);


	}
	
}


/// @brief
/// GL points and weights, first order derivative matrix
/// @param n polynomial order
/// @param k polynomial level
/// @param gl_p GL points
/// @param gl_weight GL weights
/// @param first_der first order derivative matrix
void Get_nodal_2d_storage_basis(int n, int k, double** gl_p, double** gl_weight, double** first_der){
	
	gl_p[k] = new double[n + 1];
	gl_weight[k] = new double[n + 1];
	int der_size = (n + 1) * (n + 1);
	first_der[k] = new double[der_size];
	
	// generate current gl_p, gl_w, and first_der-----------------
	GL(n, gl_p[k], gl_weight[k]);

	Mth_order_polynomial_derivative_matrix(n, 1, gl_p[k], first_der[k]);
	//--------------------------------------------------------------	
//for(int i = 0; i <= n ; ++i){
//	std::cout<< gl_p_now[i] << " " << gl_w_now[i]<<"\n";
//
//}
	// assign the address of the fist element to the input pointers
	//gl_p = &gl_p_now[0];
	//gl_weight = &gl_w_now[0];
	//first_der = &first_der_now[0];

//for(int i = 0; i <= n ; ++i){
//	std::cout<< gl_p[i] << " " << gl_weight[i]<<"\n";
//
//}
}

/// @brief
/// Lagrange interpolates on the boundaries. Alorithm 89.
/// @param n polynomial order
/// @param plevel polynomial level.
/// @param lag_l Lagrange interpolates on the left boundary.
/// @param lag_r Lagrange interpolates on the right boundary.
/// @param gl_p GL points.
void Get_nodal_2d_storage_extends(int n, int plevel, double** lag_l, double** lag_r, double** gl_p){

	double* bary = new double[n + 1];

	BARW(n, gl_p[plevel], bary);

	lag_l[plevel] = new double[n + 1];
	lag_r[plevel] = new double[n + 1];


	Lagrange_interpolating_polynomial(n, -1.0, gl_p[plevel], bary, lag_l[plevel]);
	Lagrange_interpolating_polynomial(n,  1.0, gl_p[plevel], bary, lag_r[plevel]);
	
	delete[] bary;
//std::cout<< lag_l[0][0] << "\n";
}

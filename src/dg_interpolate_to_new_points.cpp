#include "dg_interpolate_to_new_points.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_affine_map.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"
#include "dg_single_index.h"
#include <vector>
#include "dg_param.h"
#include "dg_basis.h"
#include <array>
#include <unordered_map>
#include <iostream>	// test

// forward declaration------------------------------------------------------------------------------------------
void Polynomial_interpolate_matrix(std::vector<double>& x, std::vector<double>& xi, std::vector<double>& T);

void Solutions_to_children(std::array<int, 4>& keys, int p_key);

void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start_old, int start_new, int interval);

void Form_new_set_of_points(int m, int start, std::vector<double>& y);

void Two_dir_inter(int p_key, Unit* c, std::vector<double>& T_x, std::vector<double>& T_y, int n, int m);
//--------------------------------------------------------------------------------------------------------------


/// @brief
///
/// @param keys four children key. Sequence: SW -- NW -- NE -- SE
/// @parma p_key parent's key. 
void Solutions_to_children(std::array<int, 4>& keys, int p_key){


	// pointers to the four children
	Unit* c0 = local::Hash_elem[keys[0]];
	Unit* c1 = local::Hash_elem[keys[1]];
	Unit* c2 = local::Hash_elem[keys[2]];
	Unit* c3 = local::Hash_elem[keys[3]];

	// poly order (four children are the same)
	int n = c0 -> n;
	int m = c0 -> m;

	// allocate solution space
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
		// allocate space
		c0 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));
		c1 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));
		c2 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));
		c3 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));

	}

	// four new sets of points---------------------------------------------------
	std::vector<double> xl(n + 1);	// location on parent
	std::vector<double> xr(n + 1);	// location on parent
	std::vector<double> yl(m + 1);	// location on parent
	std::vector<double> yr(m + 1);	// location on parent

	Form_new_set_of_points(n, -1.0, xl);	// left boundary is -1.0, right boundary is 0. 
	Form_new_set_of_points(n,  0.0, xr);
	Form_new_set_of_points(m, -1.0, yl);
	Form_new_set_of_points(m,  0.0, yr);
	//---------------------------------------------------------------------------

	// form interpolation matrix------------------------------------------------------
	std::vector<double> T_xl; 	// interpolation matrix
	std::vector<double> T_xr; 	// interpolation matrix
	std::vector<double> T_yl; 	// interpolation matrix
	std::vector<double> T_yr; 	// interpolation matrix

	// note: children porder == parent porder, so old points use child's poly order
	Polynomial_interpolate_matrix(nodal::gl_points[n], xl, T_xl);
	Polynomial_interpolate_matrix(nodal::gl_points[n], xr, T_xr);
	Polynomial_interpolate_matrix(nodal::gl_points[m], yl, T_yl);
	Polynomial_interpolate_matrix(nodal::gl_points[m], yr, T_yr);
	//--------------------------------------------------------------------------------

	// c0
	Two_dir_inter(p_key, c0, T_xl, T_yl, n, m);
	// c1
	Two_dir_inter(p_key, c1, T_xr, T_yl, n, m);
	// c2
	Two_dir_inter(p_key, c2, T_xr, T_yr, n, m);
	// c3
	Two_dir_inter(p_key, c3, T_xl, T_yr, n, m);
}

/// @brief
/// Use the interpolation matrix in x and y direction and obtain the interpolated solutions.
/// @param p_key parents key.
/// @param c pointer to the current child.
/// @param T_x interpolation matrix in x direction.
/// @param T_y interpolation matrix in y direction.
/// @param n polynomial order in x direction.
/// @param m polynomial order in y direction.
void Two_dir_inter(int p_key, Unit* c, std::vector<double>& T_x, std::vector<double>& T_y, int n, int m){

	std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix. 

	// y direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
		
		middle[equ] = std::vector<double> ((n + 1) * (m + 1));
		
		for(int i = 0; i <= n; ++i){

			int start = Get_single_index(i, 0, m + 1);

			// interval == 1 since we are in teh y direction
			// restriction: children inderit parent's polynomial order. 
			Interpolate_to_new_points(m + 1,  m + 1, T_y,
					local::Hash_elem[p_key] -> solution[equ], middle[equ], start, start, 1);
		}
	}

	// x direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		for(int j = 0; j <= m; ++j){

			int start = Get_single_index(0, j, m + 1);

			Interpolate_to_new_points(n + 1,  n + 1, T_x,
					middle[equ], c -> solution[equ], start, start, m + 1);
			
		}
	}


}


/// @brief
/// Map the new set of points from the reference space [-1, 1] to the space after h-refinement. 
/// @param m polynomial order.
/// @param start element left boundary coordinate.
/// @param y new set of points' coordinate. 
void Form_new_set_of_points(int m, int start, std::vector<double>& y){

	for(int j = 0; j <= m; ++j){
		
		y[j] = Affine_mapping(nodal::gl_points[m][j], start, 1.0);	// element size = 1.0
	}
}


/// @brief
/// Matrix for inerpolation between two sets of points
/// @parma x Old set of points. Size [0, n].
/// @param xi New set of points. Size [0, m]. 
/// @param T The interpolation matrix. Size [0, m] * [0, n]. Do not need to pre-allocate space. 
void Polynomial_interpolate_matrix(std::vector<double>& x, std::vector<double>& xi, std::vector<double>& T){

	int m = xi.size();
	int n = x.size();

	T = std::vector<double>(m * n);	// allocate space

	std::vector<double> w(n);	// Barycentric weights for old sets of points. Size [0, n].

	BARW(n - 1, x, w);

	for(int k = 0; k < m; ++k){

		bool row_match = false;

		for(int j = 0; j < n; ++j){

			bool same_point = Almost_equal(x[j], xi[k]);

			if(same_point){	// new point is one of the collocation points

				row_match = true;

				int index = Get_single_index(k, j, n);
				T[index] = 1.0;
			}

		}

		if(!row_match){	// if not the same point

			double s{};

			for(int j = 0; j < n; ++j){

				double t = w[j] / (xi[k] - x[j]);

				int index = Get_single_index(k, j, n);

				T[index] = t;

				s += t;
			}

			for(int j = 0; j < n; ++j){

				int index = Get_single_index(k, j, n);

				T[index] /= s;
				

			}
			

		}
	}


}

/// @brief
/// Interpolation between two sets of points by matrix multiplicaiton
/// @param m The number of new set of points.
/// @param n The number of old set of points.	
/// @param T Interpolating matrix. Size m * n.
/// @param f Solutions on the old set of points.
/// @param new_f Solusion on the new set of points. Size m. 
/// @param start_old The node number of the first value in f.
/// @param start_new The node number of the first value in new_f.
/// @param interval node index interval.
void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start_old, int start_new, int interval){

	int ni = start_new;	
	for(int i = 0; i < m; ++i){

		double t{};

		int fi = start_old;

		for(int j = 0; j < n; ++j){

			int index = Get_single_index(i, j, n);

			t += T[index] * f[fi];

//std::cout<<"f_old " << f[fi] << "\n";

			fi += interval;
		}
		
		new_f[ni] = t;


		ni += interval;

	}

//std::cout<< "\n";
}

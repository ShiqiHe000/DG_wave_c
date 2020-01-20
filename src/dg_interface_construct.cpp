#include "dg_interface_construct.h"
#include "dg_unit.h"
#include "dg_param.h"
#include "dg_single_index.h"
#include "dg_local_storage.h"
#include "dg_nodal_2d_storage.h"
#include "dg_poly_level_and_order.h"

/// @brief
/// Use Lagrange Polynomial Interplants to obtain the solution on the
/// boundaries based on the interior collocation nodes.
/// Also enforce the boundary conditions on the boundary elements  
/// This subroutine only construct interfaces in X direction. 
/// @param temp temporary pointer to each unit.
/// @param k k-th local element. 
void Construct_interface_x(Unit*& temp, int k){

	int level = Poly_order_to_level(grid::nmin, temp->n);

	for(int s = 0; s < dg_fun::num_of_equation; ++s){
		for(int j = 0; j <= (temp->m + 1); ++j){

			int position = Get_single_index_3d(0, j, s, temp->n + 1, temp->m + 1);

			int nodei = Get_single_index(j, s, dg_fun::num_of_equation);

			for(int i = 0; i <= (temp->n + 1); ++i){
	
				local::solution_int_l[k][nodei] += (nodal::lagrange_l[level][i] * (temp->solution[position]));
				local::solution_int_r[k][nodei] += (nodal::lagrange_r[level][i] * (temp->solution[position]));
				position += (temp->m + 1);
			}

		}


	}


}


// use cblas
			//int position = Get_single_index_3d(0, j, s, temp->n + 1, temp->m + 1);

			//double* temp2 = &temp->solution[position];

			//int nodei = Get_single_index(j, s, dg_fun::num_of_equation);

			//local::solution_int_l[k][nodei] = cblas_ddot(temp->n + 1, nodal::lagrange_l);

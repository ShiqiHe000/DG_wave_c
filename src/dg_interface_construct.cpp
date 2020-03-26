#include "dg_interface_construct.h"
#include "dg_unit.h"
#include <vector>
#include <unordered_map>
#include "dg_param.h"
#include "dg_basis.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"

/// @brief
/// Use Lagrange interpolants to obtain the solution on the element left and right boundaries. X direction.
/// @param temp pointer points to the current element. 
void Construct_interface_x(Unit* temp){

	int m = temp -> m;
	int n = temp -> n;

	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		temp -> solution_int_l[equ] = std::vector<double>(m + 1);
		temp -> solution_int_r[equ] = std::vector<double>(m + 1);

		for(int j = 0; j <= m; ++j ){

			std::vector<double> s_array(n + 1);

			for(int i = 0; i <= n; ++i){

				int nodei = Get_single_index(i, j, m + 1);

				s_array[i] = temp -> solution[equ][nodei];
			}

			temp -> solution_int_l[equ][j] = Interpolate_to_boundary(n, s_array, nodal::lagrange_l[n]);
			temp -> solution_int_r[equ][j] = Interpolate_to_boundary(n, s_array, nodal::lagrange_r[n]);
		
		}

	}


}

/// @brief
/// Use Lagrange interpolants to obtain the solution on the element left and right boundaries. Y direction.
/// @param temp pointer points to the current element. 
void Construct_interface_y(Unit* temp){

	int m = temp -> m;
	int n = temp -> n;

	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		temp -> solution_int_l[equ] = std::vector<double>(n + 1);
		temp -> solution_int_r[equ] = std::vector<double>(n + 1);

		for(int i = 0; i <= n; ++i ){

			std::vector<double> s_array(m + 1);

			for(int j = 0; j <= m; ++j){

				int nodei = Get_single_index(i, j, m + 1);

				s_array[j] = temp -> solution[equ][nodei];
			}

			temp -> solution_int_l[equ][i] = Interpolate_to_boundary(m, s_array, nodal::lagrange_l[m]);
			temp -> solution_int_r[equ][i] = Interpolate_to_boundary(m, s_array, nodal::lagrange_r[m]);
		
		}

	}


}

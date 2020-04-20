#include "dg_unit.h"
#include "dg_mortar_construction.h"
#include <vector>
#include "dg_param.h"
#include "dg_basis.h"
#include <cmath>	// pow
//#include <algorithm>
#include "dg_nodal_2d_storage.h"
#include "dg_interpolate_to_new_points.h"
#include "dg_single_index.h"
#include <iostream>	// test

// forward declaration ---------------------------------------------------------------------------------------
void Mortar_to_elem_interpolation(int J, int n, int level, int l_max, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar, 
				std::vector<double>& mapped_points);

void L2_projection_to_element(int J, int n, int level, int l_max, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar, 
				std::vector<double>& T);

void L2_projection_to_mortar(int J, int n, int level, int l_max, double a, double b,
			 	std::vector<double>& solution_int, std::vector<double>& psi, 
				std::vector<double>& T);
//-------------------------------------------------------------------------------------------------------------



/// @brief
/// L2 projection from element to mortar.
/// @param J Maximum polynomial order between left and right element. 
/// @parma n Element polynomial order.
/// @param level element h-refinement level.
/// @param l_max Maximum h-refinement level between two elements. 
/// @parma a Mortar coorindate offset. 
/// @parma b Mortar coorindate scaling. 
/// @parma solution_int Element interface solution. 
/// @param psi Mortar solution. 
void L2_projection_to_mortar(int J, int n, int level, int l_max, double a, double b,
			 	std::vector<double>& solution_int, std::vector<double>& psi, 
				std::vector<double>& T){

	
	if(J == n && level == l_max){	// direct copy

		// U = psi
		psi = solution_int;	
	}
	else{
		std::vector<double> mapped_points(J + 1);

		for(int i = 0; i <= J; ++i){
			mapped_points[i] = (nodal::gl_points[J][i] - a) / b;	// map GL point from mortar to element

		}

		// form interpolation matrix to interpolate value on the GL points on the element to the corresponding 
		// points facing the mortar.
		Polynomial_interpolate_matrix(nodal::gl_points[n], mapped_points, T);
		
		int start_e{};
		int start_m{};
		// interpolate the solution on GL points on the element to the points that coincode with the GL points on mortar
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

			int num{};
			for(int i = 0; i <= J; ++i){

				for(int j = 0; j <= n; ++j){
	
					psi[i + start_m] += T[num] * solution_int[start_e + j];
					++num;
				}
			}

			start_e += n + 1;
			start_m += J + 1;
		}
		

	}


}

/// @brief
/// L2 projection from element to mortar.
/// @param J Maximum polynomial order between left and right element. 
/// @parma n Element polynomial order.
/// @param level element h-refinement level.
/// @param l_max Maximum h-refinement level between two elements. 
/// @parma a Mortar coorindate offset. 
/// @parma b Mortar coorindate scaling. 
/// @parma nflux_elem Element interface numerical flux. 
/// @param nflux_mortar Mortar interface numerical flux.
/// @param mapped_points collocation points mapped from mortar to element.
void L2_projection_to_element(int J, int n, int level, int l_max, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar, 
				std::vector<double>& T){

	if(J == n && level == l_max){	// direct copy (fun and geo are conforming)

		// U = psi
		nflux_elem = nflux_mortar;	
	}
	else{	// fun or geo non-conforming

		int start_m{};
		int start_e{};

		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

			for(int i = 0; i <= n; ++i){
	
				for(int j = 0; j <= J; ++j){
				
					int nodei = Get_single_index(j , i, n + 1);
	
					nflux_elem[i + start_e] += 1.0 / b * T[nodei] * 
								(nodal::gl_weights[J][j] / nodal::gl_weights[n][i]) 
								* nflux_mortar[j + start_m];
	
				}
			}

			start_m += J + 1;
			start_e += n + 1;

		}
		
	}

}

/// @brief
/// Use Lagrange interpolation to map numerical flux from mortar to element. 
/// @param J Maximum polynomial order between left and right element. 
/// @parma n Element polynomial order.
/// @param level element h-refinement level.
/// @param l_max Maximum h-refinement level between two elements. 
/// @parma b Mortar coorindate scaling. 
/// @parma nflux_elem Element interface numerical flux. 
/// @param nflux_mortar Mortar interface numerical flux.
/// @param mapped_points collocation points mapped from mortar to element.
void Mortar_to_elem_interpolation(int J, int n, int level, int l_max, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar, 
				std::vector<double>& mapped_points){

	if(J == n && level == l_max){	// direct copy (fun and geo are conforming)

		// U = psi
		nflux_elem = nflux_mortar;	
	}
	else{	// fun or geo non-conforming

		
		std::vector<double> T;	// interpolation matrix
		Polynomial_interpolate_matrix(mapped_points, nodal::gl_points[n], T);		

		int start_m{};
		int start_e{};

		std::vector<double> middle(nflux_elem.size());
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
			Interpolate_to_new_points(n + 1, J + 1, T, nflux_mortar, middle, start_m, start_e, 1);
			start_m += J + 1;	
			start_e += n + 1;
		}

		int i{};
		for(auto& v : middle){

			
			nflux_elem[i] += v / b;

			++i;
		}
	}
}

#include "dg_unit.h"
#include "dg_mortar_construction.h"
#include <vector>
#include "dg_param.h"
#include "dg_basis.h"
#include <cmath>	// pow
#include <algorithm>
#include "dg_nodal_2d_storage.h"
#include "dg_interpolate_to_new_points.h"
#include "dg_single_index.h"
#include <iostream>	// test

// forward declaration ---------------------------------------------------------------------------------------
void Mortar_to_elem_interpolation(int J, int n, int level, int l_max, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar, 
				std::vector<double>& mapped_points);

void L2_projection_to_mortar(int J, int n, int level, int l_max, double a, double b,
			 	std::vector<double>& solution_int, std::vector<double>& psi, 
				std::vector<double>& mapped_points);

void L2_projection_to_element(int J, int n, int level, int l_max, double b,
			 	std::vector<double>& nflux_elem, std::vector<double>& nflux_mortar, 
				std::vector<double>& T);

void L2_projection_to_mortar_new(int J, int n, int level, int l_max, double a, double b,
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
				std::vector<double>& mapped_points){

	if(J == n && level == l_max){	// direct copy

		// U = psi
		psi = solution_int;	
	}
	else{

		std::vector<int> index_mortar{0, J + 1, (J + 1) * 2}; // 3 equation

		std::vector<double> bary(n + 1);

		BARW(n, nodal::gl_points[n], bary);

		mapped_points = std::vector<double>(J + 1);

//if(mpi::rank == 1){
//
//	for(auto& a : bary){
//
//		std::cout << a << " ";
//	}
//	std::cout<< "\n";
//
//}
		// L2 projection
		for(int i = 0; i <= J; ++i){	

			double s = (nodal::gl_points[J][i] - a) / b;	// map GL point from mortar to element

			mapped_points[i] = s;	// record the points (useful when mapping back form mortar to element)

			std::vector<double> lag(n + 1);

			// get the lagrange interpolation value at this point (s).
			Lagrange_interpolating_polynomial(n, s, nodal::gl_points[n], bary, lag);
//if(mpi::rank == 1){
//
////	std::cout<< "i " << i << " gl_point_mortar " << nodal::gl_points[J][i] <<" s = "<< s << "\n";
//
//	std::cout<< "s = " << s << "gl_w " << nodal::gl_weights[J][i] << "\n";
//
//	for(auto& v : lag){
//
//		std::cout << v << " ";
//	}
//	std::cout << "\n";
//}


			std::vector<int> index_elem{0, n + 1, (n + 1) * 2}; // 3 equation

			for(int j = 0; j <= n; ++j){


//if(mpi::rank == 1){
//
//	std::cout<< "j = "<< j << " gl_w "<< nodal::gl_weights[J][i] << " inter " << inter 
//				<< " lag " << lag[j] << "\n";
//}
				for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

					psi[index_mortar[equ]] += lag[j] * solution_int[index_elem[equ]];
//if(mpi::rank == 1){
//	if(equ == 1){
////		std::cout << "equ = " << equ << " psi " << psi[index_mortar[equ]] << "\n";
//
//		std::cout << "lag " << lag[j] << " gl_w " << nodal::gl_weights[J][i] << 
//				" solu " << solution_int[index_elem[equ]] << " psi " <<
//				psi[index_mortar[equ]] << "\n";
//	}
//}

				}
				std::transform(index_elem.begin(), index_elem.end(), index_elem.begin(), 
						[](int x){return x + 1;});
			}

			std::transform(index_mortar.begin(), index_mortar.end(), index_mortar.begin(), 
					[](int x){return x + 1;});
		}


	}


}


void L2_projection_to_mortar_new(int J, int n, int level, int l_max, double a, double b,
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
//for(auto& h : mapped_points){
//	std::cout<< h << "\n";
//}
		// form interpolation matrix to interpolate value on the GL points on the element to the corresponding 
		// points facing the mortar.
		Polynomial_interpolate_matrix(nodal::gl_points[n], mapped_points, T);
		
		int start_e{};
		int start_m{};
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
//for(auto& h : mapped_points){
//
//	std::cout << h << "\n";
//}
//std::cout<< "\n";
//for(auto& h : T){
//
//	std::cout<< h << "\n";
//
//}
//std::cout<< "\n";

		int start_m{};
		int start_e{};
//std::cout<< "*********************************************************************** \n";
		std::vector<double> middle(nflux_elem.size());
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
			Interpolate_to_new_points(n + 1, J + 1, T, nflux_mortar, middle, start_m, start_e, 1);
			start_m += J + 1;	
			start_e += n + 1;
		}

		int i{};
		for(auto& v : middle){

//			std::cout<< v << " ";
			
			nflux_elem[i] += v / b;
//			std::cout<< v << " " << b << "\n";

			++i;
		}
	}
}

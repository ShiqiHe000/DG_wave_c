#include <vector>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_a_times_spatial_derivative.h"
#include <unordered_map>
#include "dg_spatial_derivative.h"
#include <algorithm>
//#include "dg_param.h"


void A_times_spatial_derivative_x(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_number; ++k){

		double del_x = temp -> xcoords[1] - temp -> xcoords[0];
	
//		int n = temp -> n;
//		int m = temp -> m;
//
//		int size = (m + 1) * (n + 1);

//		std::unordered_map<int, std::vector<double>> flux_x(size);	// horizontal fluxes <equation_num, xflux>

		int index{};
		std::vector<int> index_equ{0, temp -> m + 1, (temp -> m + 1) * 2};
		for(int j = 0; j <= (temp -> m); ++j){

			std::unordered_map<int, std::vector<double>> flux_x; // horizontal fluxes <equation_num, xflux>
			std::unordered_map<int, std::vector<double>> flux_der; // <equation_num, flux_der>

			for(int i = 0; i <= (temp -> n); ++i){

				xflux(temp -> solution, flux_x, index);			
				++index;
			}
			
			// flux_der
			Spatial_derivative(temp -> n, flux_x, flux_der, temp, index_equ){

			std::transform(index_equ.begin(), index_equ.end(), index_equ.begin(),
					 [](int x){return (x + 1);});
			
		}
		


		temp = temp -> next;
	}
}


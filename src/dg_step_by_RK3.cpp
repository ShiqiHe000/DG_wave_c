#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_step_by_RK3.h"
#include "dg_time_derivative.h"
#include "dg_param.h"
#include "dg_single_index.h"

/// @brief
/// The integration in time by using low storage third order Runge-Kutta.  
/// The 3rd order Runge-Kutta methond is an explicit time time integration 
/// method. So there is a time step limitation, which depends on the
/// eigenvalues of the derivative matrix. 
/// @param tn current time
/// @param delta_t time step size
void DG_step_by_RK3(double tn, double delta_t){

	static const double am[3]{0.0, -5.0/9.0, -153.0/128.0};
	static const double bm[3]{0.0, 1.0/3.0, 3.0/4.0};
	static const double gm[3]{1.0/3.0, 15.0/16.0, 8.0/15.0};

	// equivalent to k
	
	// thrid order RK
	for(int k = 0; k < 3; ++k){

		double t = tn + bm[k] * delta_t;

		// time derivative at current time step
		DG_time_der(t);
		
		Unit* temp = local::head;
		for(int elem_k = 0; elem_k < local::local_elem_num; ++elem_k){

			std::unordered_map<int, std::vector<double>> G;

			for(int l = 0; l < dg_fun::num_of_equation; ++l){

				int size = (temp -> n + 1) * (temp -> m + 1);

				G[l] = std::vector<double>(size);

				for(int j = 0; j <= (temp -> m); ++j){

					for(int i = 0; i <= (temp -> n); ++i){

						int nodei = Get_single_index(i, j, temp -> m + 1);

						G[l][nodei] = am[k] * G[l][nodei] + (temp -> solution_time_der)[l][nodei];

						(temp -> solution)[l][nodei] += gm[k] * delta_t * G[l][nodei];
	
					}

				}
			

			}
			(temp -> solution_time_der).clear();
			temp = temp -> next;
		}

	}
}

#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_step_by_RK3.h"

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
		
		
		Unit* temp = local::head;

	}
}

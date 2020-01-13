#include "dg_end_game.h"
#include "dg_poly_level_and_order.h"
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"

/// @brief
/// Free memory on the heap
void Free_local(){

	

	int level_max = Poly_order_to_level(grid::nmin, grid::nmax);

	for(int k = 0; k <= level_max; ++k){

		delete[] nodal::gl_p[k];
		delete[] nodal::gl_w[k];
		delete[] nodal::first_der[k];
		delete[] nodal::lagrange_l[k];
		delete[] nodal::lagrange_r[k];

	}	

	delete[] nodal::gl_p;
	delete[] nodal::gl_w;
	delete[] nodal::first_der;
	delete[] nodal::lagrange_l;
	delete[] nodal::lagrange_r;

}

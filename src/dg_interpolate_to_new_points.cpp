#include "dg_interpolate_to_new_points.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_affine_map.h"
#include "dg_nodal_2d_storage.h"

/// @brief
///
/// @param key child's key
void Solutions_on_child(int key){

	Unit* temp = local::Hash_elem[key];	// pointer to this elem

	double del_x = (temp -> xcoords[1]) - (temp -> coords[0]);

	for(int j = 0; j <= (temp -> m); ++j){

		std::vector<double> x(temp -> n + 1);	// location on parent

		// generate new sets of poitns
		for(int i = 0; i <= (temp -> n); ++i){

			x[i] = Affine_mapping(nodal::gl_points[temp -> n][i], temp -> xcoords[0], del_x);
			
		}

		// interpolate to new sets of points

	}

}

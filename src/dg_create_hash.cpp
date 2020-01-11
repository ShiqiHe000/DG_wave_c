#include "dg_create_hash.h"
#include "dg_local_storage.h"
#include <unordered_map>
#include "dg_unit.h"
#include "dg_single_index.h"
#include "dg_hilbert_curve.h"
#include "dg_cantor_pairing.h"
#include "dg_param.h"
#include "dg_index_local_global.h"

/// @brief
/// Create the hash table <key, unit_elem>.
/// key is generated by the element index(i, j, k).
/// unit_elem is the element unit storage.
void Create_hash(){
	
	// create hash table
//	std::unordered_map<int, Unit*> Hash_elem;

	int key_pre{};
	int solution_num = dg_fun::num_of_equation * (grid::nmin + 1) * (grid::nmin + 1);

	// build each unit
	for(int k = 0; k < local::local_elem_num; ++k){

		// 2d integer coordinates
		int g_index = Index_local_to_global(mpi::rank, k);
		int ii, jj;
		d2xy(grid::exp_x, g_index, jj, ii );

		int key = Get_key_fun(ii, jj, 0);
		
		// create a unit
		local::Hash_elem[key] = new Unit();
		
		// if first element, set head ptr
		if(k == 0){

			local::head = local::Hash_elem[key];	//head points to the first unit 
		}		
		
		// index
		local::Hash_elem[key] -> index[0] = ii;
		local::Hash_elem[key] -> index[1] = jj;

		// status
		local::Hash_elem[key] -> status = local::status[k];
		
		// coordinates
		for(int i = 0; i < 2; ++i){
			int nodei = Get_single_index(k, i, 2);
			local::Hash_elem[key] -> xcoords[i] = local::x_local[nodei];
			local::Hash_elem[key] -> ycoords[i] = local::y_local[nodei];
		}


		// solution
		local::Hash_elem[key] -> solution = new double[solution_num]{};
	
		// previous unit should point to current unit
		if(k > 0){
			local::Hash_elem[key_pre] -> next = local::Hash_elem[key];
		}

		key_pre = key;
	}

	// free memory
	delete[] local::x_local;
	delete[] local::y_local;
}
#include "dg_create_hash.h"
#include "dg_local_storage.h"
#include <unordered_map>
#include "dg_unit.h"
#include "dg_single_index.h"
#include "dg_hilbert_curve.h"
#include "dg_cantor_pairing.h"
#include "dg_param.h"
#include "dg_index_local_global.h"
#include <iostream>	// test

/// @brief
/// Create the hash table <key, unit_elem>.
/// key is generated by the element index(i, j, k).
/// unit_elem is the element unit storage.
void Create_hash(){
	
	// create hash table
	std::unordered_map<int, Unit*> Hash_elem;

	int key_pre{};

	// build each unit
	for(int k = 0; k < local::local_elem_num; ++k){

		// 2d integer coordinates
		int g_index = Index_local_to_global(mpi::rank, k);
		int ii, jj;
		d2xy(grid::exp_x, g_index, jj, ii );

		int key = Get_key_fun(ii, jj, 0);
//std::cout<< mpi::rank<< " " << ii <<" "<< jj<< " " << key << "\n";
		// create a unit
		Hash_elem[key] = new Unit();
//		std::cout<< Hash_elem[k]->level <<" " <<Hash_elem[k] -> n<< "\n";
		
		// if first element, set head ptr
		if(k == 0){

			local::head = Hash_elem[key];	//head points to the first unit 
		}		
		
		// index
		Hash_elem[key] -> index[0] = ii;
		Hash_elem[key] -> index[1] = jj;

		// status
		Hash_elem[key] -> status = local::status[k];
//std::cout << mpi::rank << " " << k << " " << Hash_elem[key]->index[1] <<" " << Hash_elem[key] -> status << "\n";
		// coordinates
		for(int i = 0; i < 2; ++i){
			int nodei = Get_single_index(k, i, 2);
			Hash_elem[key] -> xcoords[i] = local::x_local[nodei];
			Hash_elem[key] -> ycoords[i] = local::y_local[nodei];
//			std::cout<< mpi::rank << " "<<Hash_elem[k]->xcoords[i] << " " << Hash_elem[k]->ycoords[i] << "\n";
		}
	
		// previous unit should point to current unit
		if(k > 0){
			Hash_elem[key_pre] -> next = Hash_elem[key];
		}

		key_pre = key;
	}

	// free memory
	delete[] local::x_local;
	delete[] local::y_local;
}

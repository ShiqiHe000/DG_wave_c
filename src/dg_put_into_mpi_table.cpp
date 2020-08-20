#include <vector>
#include <unordered_map>
#include "dg_unit.h"
#include "dg_put_into_mpi_table.h"
#include "dg_boundary_table.h"
#include "dg_local_storage.h"
#include "dg_cantor_pairing.h"
#include <algorithm>

/// @brief 
/// Insert the current MPI boundary element to the MPI boundary table. 
/// @param temp Pointer to the current unit. 
/// @param facen_it Iterator on the element face info vector.
/// @param mpi_table The relevent MPI bountary table. 
/// @param target_rank The neighbour's rank. 
void Put_in_mpi_table(Unit* temp, std::vector<Unit::Face>::iterator& facen_it, 
			std::unordered_map<int, std::vector<mpi_table>>& table){

	if(table.count(facen_it -> rank) == 0){	// if this rank has not been not record yet

		table[facen_it -> rank] = std::vector<mpi_table>();

		long long int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		// mpi_length and owners_rank will be recorded later
		table[facen_it -> rank].push_back({local_key, 0, 0});
	}
	else{ // this rank is already been record

		long long int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		// avoid duplication
		auto it = std::find_if(table[facen_it -> rank].begin(), table[facen_it -> rank].end(),
					[local_key] (const mpi_table& v) {return v.local_key == local_key;});

		if(it == table[facen_it -> rank].end()){	// if not find

			table[facen_it -> rank].push_back({local_key, 0, 0});
		}

	}

}


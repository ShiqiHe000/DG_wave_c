#include <iostream>
#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"
#include "dg_write_mpi_table.h"

void Write_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& table){

	for(auto& v : table){

		int target_rank = v.first;

		std::cout << "target_rank " << target_rank << "\n";

		for(auto it = v.second.begin(); it != v.second.end(); ++it){

			std::cout<< "local_key " << it -> local_key << " m_length " << it -> mpi_length << "\n";

		}

	}
	std::cout << "===========================\n";
		
}


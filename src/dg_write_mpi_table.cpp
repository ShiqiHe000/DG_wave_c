#include <iostream>
#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include "dg_write_mpi_table.h"
#include "dg_param.h"
#include <mpi.h>

int file_table = 1;
//-------------------------------------------------------------------------------
void Write_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& table1, 
			std::unordered_map<int, std::vector<mpi_table>>& table2);

void Write_table_all(std::unordered_map<int, std::vector<mpi_table>>& table1, 
			std::unordered_map<int, std::vector<mpi_table>>& table2);
// ------------------------------------------------------------------------------

void Write_table_all(std::unordered_map<int, std::vector<mpi_table>>& table1, 
			std::unordered_map<int, std::vector<mpi_table>>& table2){

	for(int i = 0; i < mpi::num_proc; ++i){

		if(mpi::rank == i){
			Write_mpi_table(table1, table2);

		}
		else{

			MPI_Barrier(MPI_COMM_WORLD);
		}

	}

}

void Write_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& table1, 
			std::unordered_map<int, std::vector<mpi_table>>& table2){


		// generate the file name
	std::stringstream ss;
	ss << "../table_owner/table" << std::setfill('0') << std::setw(5) << file_table << ".dat";
	std::string filename = 	ss.str();
	std::ofstream myfile; 	// stream class to write on files	

	if(mpi::rank == 0){
		myfile.open(filename, std::ios::trunc); // truncate the old file
	}
	else{
		myfile.open(filename, std::ios::app);

	}

	myfile<< "===============" << mpi::rank << "==============================="<< "\n";
	myfile << "*********** table 1 *************************************** \n";

	int tol_length{};

	for(auto& v : table1){

		tol_length = 0;

		int target_rank = v.first;
					
		myfile<<"target_rank " << target_rank << " ----------------------------------" << "\n";
	
		for(auto it = v.second.begin(); it != v.second.end(); ++it){

			myfile << "local_key " << it -> local_key << " m_length " << it -> mpi_length <<
				" owner " << it -> owners_rank << "\n";

			tol_length += it -> mpi_length;
		}
				
		myfile<< "tol_length = "<< tol_length << "\n";
	}	


	myfile << "*********** table 2 *************************************** \n";
	for(auto& v : table2){

		tol_length = 0;

		int target_rank = v.first;
					
		myfile<<"target_rank " << target_rank << " ----------------------------------" << "\n";
	
		for(auto it = v.second.begin(); it != v.second.end(); ++it){

			myfile << "local_key " << it -> local_key << " m_length " << it -> mpi_length <<
				" owner " << it -> owners_rank << "\n";

			tol_length += it -> mpi_length;
		}
				
		myfile<< "tol_length = "<< tol_length << "\n";

	}	

	myfile.close();
	file_table++;
	
}


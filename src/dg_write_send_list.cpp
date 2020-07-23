#include <iostream>
#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include "dg_write_send_list.h"
#include "dg_param.h"
#include <mpi.h>
#include "dg_local_storage.h"

int file_send = 1;

void Write_send_list();
void Write_send_list_all();

void Write_send_list_all(){

	for(int i = 0; i < mpi::num_proc; ++i){

		if(mpi::rank == i){
			Write_send_list();

		}
		else{

			MPI_Barrier(MPI_COMM_WORLD);
		}

	}


}

void Write_send_list(){

		// generate the file name
	std::stringstream ss;
	ss << "../send_list/sendlist" << std::setfill('0') << std::setw(5) << file_send << ".dat";
	std::string filename = 	ss.str();
	std::ofstream myfile; 	// stream class to write on files	

	if(mpi::rank == 0){
		myfile.open(filename, std::ios::trunc); // truncate the old file
	}
	else{
		myfile.open(filename, std::ios::app);

	}

	myfile<< "=============== my rank " << mpi::rank << "==============================="<< "\n";
	myfile << "*********** pre rank "<< mpi::rank - 1 << " ****************** \n" ;


	for(auto& v : LB::Send.pre){

		myfile << v << "\n";
				
	}	


	myfile << "*********** next rank "<< mpi::rank + 1 << " ****************** \n" ;
	for(auto& v : LB::Send.next){

		myfile << v << "\n";
				
	}	

	myfile.close();
	file_send++;	
}

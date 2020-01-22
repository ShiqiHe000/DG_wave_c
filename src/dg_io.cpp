#include "dg_io.h"
#include <mpi.h>
#include "dg_param.h"
#include "dg_local_storage.h"
#include <iomanip>	// std::setfill
#include <string>
#include <iostream>
#include <sstream>	// std::stringstream
#include <fstream>	// read and write to file
#include "dg_unit.h"


// forward declaration
void Write_mesh();

/// @brief
/// Output data in serial order.
/// @param t current time.
void Serial_io(double t){

	for(int k = 0; k < mpi::num_proc; ++k){

		if(mpi::rank == k){
			Write_mesh();

		}
		else{

			MPI_Barrier(MPI_COMM_WORLD);		
		
		}
		

	}
	

}

// global file number
int file_num = 1;


/// @brief
/// Processor 1 create the file. All the processor write the data in order.
/// In the end processor 1 close the file.
/// 
void Write_mesh(){

	// generate the file name
	std::stringstream ss;
	ss << "../outputs/output" << std::setfill('0') << std::setw(5) << file_num << ".dat";
	std::string filename = 	ss.str();
	std::ofstream myfile; 	// stream class to write on files	

	// traverse the linked-list
	Unit* temp = local::head;

	// processor open file
	if(mpi::rank == 0){
		
		myfile.open(filename, std::ios::trunc);	// truncate the old file

		// headers
		myfile<< "TITLE = \"MESH AND SOLUTIONS\" \n";
		myfile<< "VARIABLES = \"X\", \"Y\", \"RANK\" \n";

	}
	else{

		myfile.open(filename, std::ios::app);	// All output operations are performed at the end of the file
	}

	// write solutions
	int elem = 1;

	for(int iel = 0; iel < local::local_elem_num; ++iel){

		myfile << "ZONE T= " << "\"" << "IEL" << std::setw(6) << elem << "\"" << "  " 
			<<"I=2, J=2"<< "  " << "DATAPACKING = POINT"<< "\n";
		
		++elem;
		
		myfile << std::fixed;
		myfile << std::setprecision(5);
		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[0] << "  " << mpi::rank << "\n";
		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[1] << "  " << mpi::rank << "\n";
		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[0] << "  " << mpi::rank << "\n";
		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[1] << "  " << mpi::rank << "\n";
	
		temp = temp -> next;
	}


	// close the file
	myfile.close();
	++file_num;

}

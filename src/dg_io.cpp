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
void Write_mesh(double t, int pre_elem);

/// @brief
/// Output data in serial order.
/// @param t current time.
void Serial_io(double t){

	int pre_elem{};

	MPI_Exscan(&local::local_elem_num, &pre_elem, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	for(int k = 0; k < mpi::num_proc; ++k){

		if(mpi::rank == k){
			Write_mesh(t, pre_elem);

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
/// @param t current time step
void Write_mesh(double t, int pre_elem){

	// generate the file name
	std::stringstream ss;
	ss << "../outputs/output" << std::setfill('0') << std::setw(5) << file_num << ".dat";
	std::string filename = 	ss.str();
	std::ofstream myfile; 	// stream class to write on files	

	// traverse the linked-list
	Unit* temp = local::head;
	int elem = pre_elem;

	// processor open file
	if(mpi::rank == 0){
		
		myfile.open(filename, std::ios::trunc);	// truncate the old file

		// headers
		myfile<< "TITLE = \"MESH AND SOLUTIONS\" \n";
		myfile<< "VARIABLES = \"X\", \"Y\", \"RANK\", \"HLEVEL\", \"VAR\", \n";
		
	}
	else{

		myfile.open(filename, std::ios::app);	// All output operations are performed at the end of the file
	}

	// write solutions
	for(int iel = 0; iel < local::local_elem_num; ++iel){

		myfile << std::fixed;
		myfile << std::setprecision(5);
//		myfile << "ZONE T= " << "\"" << "IEL" << std::setw(6) << elem << "\"" << "  " 
//			<<"I=2, J=2"<< "  " << "DATAPACKING = POINT"<< "\n";
		myfile << "ZONE T= " << "\"" << "IEL" << std::setw(6) << elem << "\"," << "  " 
			<<"I=2, J=2, "<< "SOLUTIONTIME=" << t <<", DATAPACKING = POINT" << "\n";
		
		++elem;
		
	//	myfile << std::fixed;
	//	myfile << std::setprecision(5);
		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[0] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  " <<temp -> var <<"\n";
		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[1] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  "<< temp -> var<<"\n";
		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[0] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  "<< temp -> var<< "\n";
		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[1] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  "<< temp -> var<<"\n";
		
		temp = temp -> next;
	}


	// close the file
	myfile.close();
	++file_num;

}

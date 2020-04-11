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
#include "dg_cantor_pairing.h"
#include <vector>
#include <unordered_map>
#include "dg_interface_construct.h"
#include "dg_nodal_2d_storage.h"
#include "dg_basis.h"

// forward declaration-----------------------------------------------------------------------------------
void Write_mesh(double t, int pre_elem);

void Interpolate_to_four_corner(Unit* temp, std::unordered_map<int, std::vector<double>>& four);
//-------------------------------------------------------------------------------------------------------

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
		myfile<< "VARIABLES = \"X\", \"Y\", \"RANK\", \"HLEVEL\", \"KEY\", \"PRESSURE\", \"U\", \"V\",\n";
		
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
		
		// interpolate solution to the four corner----------------------------------------------------
		std::unordered_map<int, std::vector<double>> four;
		Interpolate_to_four_corner(temp, four);
		//--------------------------------------------------------------------------------------------


	//	myfile << std::fixed;
	//	myfile << std::setprecision(5);
		int key_now = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[0] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  "<< key_now 
			<< "  " << four[0][0]<< "  " << four[1][0] << "  "<< four[2][0] <<"\n";


		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[1] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  "<< key_now
			<< "  " << four[0][1] <<"  " << four[1][1] <<"  "<< four[2][1]<<"\n";


		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[0] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  "<<key_now 
			<< "  " << four[0][2] << "  "<< four[1][2] << "  "<< four[2][2]<<"\n";


		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[1] 
			<< "  " << mpi::rank << "  " << temp -> index[2]<< "  "<<key_now
			<< "  " << four[0][3] << "  " << four[1][3] << "  " << four[2][3]<<"\n";
		
		temp = temp -> next;
	}


	// close the file
	myfile.close();
	++file_num;

}

/// @brief
/// Interpolates the interior solutions to four corner points
/// @param temp pointer to the current element.
/// @param four solutions on the four corner points. <equ, four solutions>
// four solution sequnece
//     2 --------- 3
//     |           |
//     |           |
//     |           |
//     |           |
//     0 --------- 1
void Interpolate_to_four_corner(Unit* temp, std::unordered_map<int, std::vector<double>>& four){

	// construct interface in y direction.
	Construct_interface_y(temp);

	int a{};
	
	// then interpolates to the four corners
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		// copy the vector
		std::vector<double> s_left(temp -> n + 1);
		std::vector<double> s_right(temp -> n + 1);

		for(int i = 0; i <= (temp -> n); ++i){

			s_left[i] = temp -> solution_int_l[a];
			s_right[i] = temp -> solution_int_r[a];
			++a;	

		}

		// allocate space for corner hash
		four[equ] = std::vector<double> (4);

		// lower left
		four[equ][0] = Interpolate_to_boundary(temp -> n, s_left, nodal::lagrange_l[temp -> n]);		

		// lower right
		four[equ][1] = Interpolate_to_boundary(temp -> n, s_right, nodal::lagrange_l[temp -> n]);		

		// upper left
		four[equ][2] = Interpolate_to_boundary(temp -> n, s_left, nodal::lagrange_r[temp -> n]);		

		// upper right
		four[equ][3] = Interpolate_to_boundary(temp -> n, s_right, nodal::lagrange_r[temp -> n]);		

	}


	// clear solution_int_l/r
	(temp -> solution_int_l).clear();
	(temp -> solution_int_r).clear();
}



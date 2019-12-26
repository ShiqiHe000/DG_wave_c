#include <mpi.h>
#include <iostream>
#include <fstream> 	// I/O with files
#include <string>
#include <sstream> 	// mainipulate a string
#include "dg_param.h"	
#include "dg_read_mesh_2d.h"

/// convert the 2d array index into a single 1d array index
/// @param row row number (start with 0) 
/// @param col column number (start with 0)
/// @param total_col total number of column (start with 1)
int Get_single_index(int row, int col, int total_col){

	return (row * total_col + col);

}

void Read_mesh_2d(){

	if(mpi::rank == 0){
	
		std::cout<< "------------------------------" << "\n";
		std::cout<< "Start to read mesh." << "\n";
		std::cout<< "------------------------------" << "\n";
		
		// read from file
		std::ifstream gmsh_file;
		
		// open the file
		gmsh_file.open(fileinfo::fileplace);

		if(gmsh_file.is_open()){
			std::cout << "file open"<< "\n";
		}

		// store the information in a line here
		std::string charline;
	
		// skip dummy lines
		const int to_node = 12;
		for(int i = 0; i < 12; i++ ){
			
			std::getline(gmsh_file, charline);
		
		
		}	
		
		// now we are at "$Nodes"
		std::stringstream sss(charline);
		std::string string1;
		sss >> string1;

		// check if we are at the target position
		if(string1 != "$Nodes"){
			std::cout<< "Problem in .msh file. Please use Legency ASCII2 format." << "\n";
		}
		
		// get the total node number
		std::getline(gmsh_file, charline);
		
		std::stringstream ss(charline);
		int total_node;
		ss >> total_node;
		
		// allocate a memory space on heap to hold the 2d array
		const int node_size = 2 * total_node; 
		double* ptr_node;	// pointer to node coordinate
		ptr_node = new double[node_size]; 	// store the node coordinates
	
		

//		int a;
//		for(int i = 0; i < num_of_phy_name; i++){
//			
//			std::getline(gmsh_file, charline);
//			std::stringstream ss(charline);
//			ss >> a;
//		}
//
//		if(a != 2){
//			std::cout<< "You forgot to set the physical surface." << "\n";
//			std::cout<< charline << "\n";
//			return;
//		}



		// close file
		gmsh_file.close();
	
	
	
	}

}


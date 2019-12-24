#include <mpi.h>
#include <iostream>
#include <fstream> 	// I/O with files
#include <string>
#include <sstream> 	// mainipulate a string
#include "dg_param.h"	
#include "dg_read_mesh_2d.h"

// Trim a string
std::string Rtrim(std::string& string1);

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

		std::cout<< string1 << "\n";
		if(string1 != "$Nodes"){
			std::cout<< "Problem in .msh file. Please use Legency ASCII2 format." << "\n";
		}
		
		// get the total node number
		std::getline(gmsh_file, charline);
		
		std::stringstream ss(charline);
		int total_node;
		ss >> total_node;
		
		double* ptr_node;	// pointer to node coordinate
		prt_node = new double[2, total_node];
		
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

// define right trim a string function
std::string Rtrim(std::string& string1){
	
	const std::string targets = " ";

	string1.erase(string1.find_last_not_of(targets)+1);

	return string1;

}

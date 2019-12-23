#include <mpi.h>
#include <iostream>
#include <fstream> 	// I/O with files
#include <string>
#include <sstream> 	// mainipulate a string
#include "dg_param.h"	
#include "dg_read_mesh_2d.h"

void Read_mesh_2d(){

	if(mpi::rank == 0){
	
		std::cout<< "------------------------------" << "\n";
		std::cout<< "Start to read mesh." << "\n";
		std::cout<< "------------------------------" << "\n";
		
		// read from file
		//std::string* path = fileinfo::fileplace;
		std::ifstream gmsh_file();
		
		// open the file
		gmsh_file.open("../gmsh_files/example.txt");

		if(gmsh_file.is_open()){
			std::cout << "file open"<< "\n";
		}

		// store the information in a line here
		std::string charline;

		while(std::getline(gmsh_file, charline)){	// while there is a line
			std::istringstream ss(charline);

			int a, b, c;

			ss >> a >> b >> c;

			std::cout<< a << " " << b << " " << c << "\n";	
		
		}



		// close file
		gmsh_file.close();
	
	
	
	}

}


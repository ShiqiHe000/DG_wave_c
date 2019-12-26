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

///
void Read_mesh_2d(){

	if(mpi::rank == 0){
	
		std::cout<< "------------------------------" << "\n";
		std::cout<< "Start to read mesh." << "\n";
		std::cout<< "------------------------------" << "\n";
		
		// read from file
		std::ifstream gmsh_file;
		
		// open the file
		gmsh_file.open(fileinfo::fileplace);

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
		double* node_xy;	// pointer to node coordinate
		node_xy = new double[node_size]; 	// store the node coordinates

		// read the node coordinates
		for(int i = 0; i < total_node; ++i){
			std::getline(gmsh_file, charline);
			std::stringstream coord(charline);
			
			int nodei;

			coord >> nodei >> node_xy[Get_single_index(i, 0, 2)]
				>> node_xy[Get_single_index(i, 1, 2)];
		
		}	

		// read dummy lines
		std::getline(gmsh_file, charline);	// $EndNodes
		std::getline(gmsh_file, charline);	// $Elements

		// number of elements
		std::getline(gmsh_file, charline);
		std::stringstream elem(charline);
		int num_of_element;
		elem >> num_of_element;

		// record 4 nodes number of each element
		const int quad_size = 4 * num_of_element;
		double* quad_node = new double[quad_size];
		int total_quad{0};	// total quard-element number
		
		for(int i = 0; i < num_of_element; ++i){
			std::getline(gmsh_file, charline);

			int a, elem_type;	// get element type

			std::stringstream type(charline);

			type >> a >> elem_type;

			if(elem_type == 3){

				int b, c;

				type >> a >> b >> c  
					>> quad_node[Get_single_index(total_quad, 0, 4)] 
					>> quad_node[Get_single_index(total_quad, 1, 4)]
					>> quad_node[Get_single_index(total_quad, 2, 4)]
					>> quad_node[Get_single_index(total_quad, 3, 4)];	
				
				++total_quad;
			}

		
		}

		
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


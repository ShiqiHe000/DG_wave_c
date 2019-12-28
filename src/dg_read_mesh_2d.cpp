#include <mpi.h>
#include <iostream>
#include <fstream> 	// I/O with files
#include <string>
#include <sstream> 	// mainipulate a string
#include <cmath>	// std::abs()
#include "dg_param.h"	
#include "dg_read_mesh_2d.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"

// forward declaration----------------------------------------------------------------------------
void Sort_node_ordering(int &total_node, int &total_quad, int*& quad_node, double*& node_xy);


void Get_standard(double &x1, double &y1, 
		double &x2, double &y2,
		double &x_max, double &x_min, 
		double &y_max, double y_min);

int Get_scores(double& x_max, double& y_max, double& node_xy1, double& node_xy2, 
		const int (&x_scores)[2], const int (&y_scores)[2]);

bool AlmostEqual(double& a, double& b);
//-------------------------------------------------------------------------------------------------


/// @brief
/// Read .msh file, record the xy coordinates of each elements. 
/// Sorting the node-ordering format. Warrent each element is ordering 
///  as: 
/// \verbatim
///  4--------3  
///  | 	      | 
///  |	      | 
///  |	      |  
///  2--------1
/// \endverbatim
/// Get the ij coordinate.\n
/// @note the .msh file has to be in version 2 ASII format. 
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
		int* quad_node = new int[quad_size];
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
		
		// Sort the node ordering format to the standard way
		Sort_node_ordering(total_node, total_quad, quad_node, node_xy );


		// free memory on heap
		delete[] quad_node;
		delete[] node_xy;

		// close file
		gmsh_file.close();
	
	}

}


/// Sprting the node-ordering format. Warrent each element is ordered as below: \n
/// \verbatim
///  4--------3  
///  | 	      | 
///  |	      | 
///  |	      |  
///  2--------1
/// \endverbatim
void Sort_node_ordering(int &total_node, int &total_quad, int*& quad_node, double*& node_xy ){
	
	//-----------------------------------
	const int x_scores[]{2, 1};
	const int y_scores[]{20, 10};

	const int node1 = x_scores[1] + y_scores[1];
	const int node2 = x_scores[1] + y_scores[0];
	const int node3 = x_scores[0] + y_scores[0];
	const int node4 = x_scores[0] + y_scores[1];
	//----------------------------------
	
	//------------------------------------------
	SortMesh::elem_x_position = new double[4 * total_quad]{}; 
	SortMesh::elem_y_position = new double[4 * total_quad]{}; 
	//------------------------------------------
	

	for(int k = 0; k < total_quad; ++k){

		double x_max, x_min, y_max, y_min;
		
		int node11, node33;

		node11 = quad_node[Get_single_index(k, 0, 4)] - 1;	// node id start with 0
		node33 = quad_node[Get_single_index(k, 2, 4)] - 1;
		
		Get_standard(node_xy[Get_single_index(node11, 0, 2)],
			       node_xy[Get_single_index(node11, 1, 2)], 
			       node_xy[Get_single_index(node33, 0, 2)], 
			       node_xy[Get_single_index(node33, 1, 2)], 
		       		x_max, x_min, y_max, y_min);	

		for(int i = 0; i < 4; ++i){
			
			int inode = quad_node[Get_single_index(k, i, 4)] - 1;
			
			int score = Get_scores(x_max, y_max, node_xy[Get_single_index(inode, 0, 2)], 
					node_xy[Get_single_index(inode, 1, 2)], 
					x_scores, y_scores);
			

			if(score == node1){

				int ii = Get_single_index(k, 0, 4);
				SortMesh::elem_x_position[ii] = node_xy[Get_single_index(inode, 0, 2)];
				SortMesh::elem_y_position[ii] = node_xy[Get_single_index(inode, 1, 2)];	

			}
			else if(score == node2){
			
				int ii = Get_single_index(k, 1, 4);
				SortMesh::elem_x_position[ii] = node_xy[Get_single_index(inode, 0, 2)];
				SortMesh::elem_y_position[ii] = node_xy[Get_single_index(inode, 1, 2)];	
			}
			else if(score == node3){
				
				int ii = Get_single_index(k, 2, 4);
				SortMesh::elem_x_position[ii] = node_xy[Get_single_index(inode, 0, 2)];
				SortMesh::elem_y_position[ii] = node_xy[Get_single_index(inode, 1, 2)];	
			}
			else if(score == node4){
				
				int ii = Get_single_index(k, 3, 4);
				SortMesh::elem_x_position[ii] = node_xy[Get_single_index(inode, 0, 2)];
				SortMesh::elem_y_position[ii] = node_xy[Get_single_index(inode, 1, 2)];	
			
			} 
			else{
				std::cout<< "In dg_read_2d_mesh.cpp, the element 'score' doesn not \
					match the stardards." << "\n";
				return;
			}
		}

	}

	
}

void Get_standard(double &x1, double &y1, 
		double &x3, double &y3,
		double &x_max, double &x_min, 
		double &y_max, double y_min){

	if(x1 > x3){
		x_max = x1;
		x_min = x3;
	}
	else{
		x_max = x3;
		x_min = x1;
	
	}
	
	if(y1 > y3){
		y_max = y1;
		y_min = y3;
	}
	else{
		y_max = y3;
		y_min = y1;
	
	}

}

int Get_scores(double& x_max, double& y_max, double& node_xy1, double& node_xy2, 
		const int (&x_scores)[2], const int (&y_scores)[2]){
	
	int score = 0;
	bool flag1, flag2;

	flag1 = AlmostEqual(node_xy1, x_max);
	flag2 = AlmostEqual(node_xy2, y_max);

	if(flag1){
		score += x_scores[0];
	
	}
	else{
	
		score += x_scores[1];
	}
	
	if(flag2){
		score += y_scores[0];
	}
	else{
		score += y_scores[1];
	}
	
	return score;
	

}


/// @brief
/// If the two double precision number a and b are equal (|a - b| <= 1.0e-11),
/// then return true.
bool AlmostEqual(double& a, double& b){
	
	const double threshold = 1.0e-11;

	if(std::abs(a - b) <= threshold){
		return true;
	}
	else{
		return false;
	}

}

#include <iostream>
#include <mpi.h>
#include "dg_param.h"
#include "dg_cross_section.h"
#include <cassert>
#include <sstream>	// std::stringstream
#include <fstream>	// read and write to file
#include <string>
#include "dg_basis.h"
#include "dg_affine_map.h"
#include "dg_single_index.h"
#include <vector>
#include "dg_unit.h"
#include "dg_nodal_2d_storage.h"
#include <unordered_map>
#include "dg_local_storage.h"

// forward declaration -------------------------------------
void Solution_cross_section(double x);
void Write_cross_section(double x_fixed);
// ---------------------------------------------------------


/// @brief
/// Outout the solution on the cross-seciton of the wave. 
/// Note the cross section is vertical to the x axis. 
void Solution_cross_section(double x){

	assert(x <= grid::gx_r && x >= 0.0 && "The cross-section is out of the bound.\n");

	for(int k = 0; k < mpi::num_proc; ++k){

		if(mpi::rank == k){

			Write_cross_section(x);

		}
		else{
			MPI_Barrier(MPI_COMM_WORLD);
		}
	

	}
	
}

/// @brief
/// Interpolate the solutions to the cross-section and write to the file. 
/// @param x_fixed x coordinate (cross-section is perpendicular to x axis). 
/// @param elem element number on the cross section. 
void Write_cross_section(double x_fixed){

	// generate the file name
	std::string filename = fileinfo::crosssection_filename;
	std::ofstream myfile; 	// stream class to write on files	

	// traverse the linked-list
	Unit* temp = local::head;

	// processor open file
	if(mpi::rank == 0){
		
		myfile.open(filename, std::ios::trunc);	// truncate the old file

		// headers
		myfile<< "y pressure \n";
		
	}
	else{

		myfile.open(filename, std::ios::app);	// All output operations are performed at the end of the file
	}

	// write solutions
	for(int iel = 0; iel < local::local_elem_num; ++iel){

		// cross-section passes this element
		if(temp -> xcoords[0] <= x_fixed && temp -> xcoords[1] > x_fixed){
			
			double del_x = temp -> xcoords[1] - temp -> xcoords[0];
			double del_y = temp -> ycoords[1] - temp -> ycoords[0];
			
			// coordinate on reference sapce (x direction)
			double x_ref = Map_to_reference(x_fixed, temp -> xcoords[0], del_x);	

			// barycentric weights
			std::vector<double> bary(temp -> n + 1);
			BARW(temp -> n, nodal::gl_points[temp -> n], bary);

			// interpolate pressure to this cross-section
			for(int j = 0; j <= temp -> m; ++j){
					
				std::vector<int> index(temp -> n + 1);
				
				for(int i = 0; i < temp -> n + 1; ++i){
	
					index[i] = Get_single_index(i, j, temp -> m + 1);
				}
				double res = Lagrange_interpolation(temp -> n, x_ref, nodal::gl_points[temp -> n], 
						temp -> solution[0], bary, index);

				// y coordinate, pressure
				double y = Affine_mapping(nodal::gl_points[temp -> m][j], temp -> ycoords[0], del_y);
				myfile << y << " " << res << "\n";
	
			}

		}
	
		temp = temp -> next;
	
	}	

	myfile.close();

}

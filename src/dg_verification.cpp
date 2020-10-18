#include "dg_local_storage.h"
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_unit.h"
#include "dg_user_defined.h"
#include <cmath>	// std::abs
#include "dg_verification.h"
#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>	// std::stringstream
#include <fstream>	// read and write to file
#include <iomanip>	// std::setfill
#include <iostream>	// test

// global elem -----------------
int g_elem{};
// ----------------------------

// ----------------------------------------------
void Get_error();
void Write_error();
void Write_error_each_proc(double t);
// ----------------------------------------------

/// @brief
/// Verfies your results. Now assume uniform grids. 
void Get_error(){

	result::L2_norm = std::vector<double>(dg_fun::num_of_equation);

	// temperary pointer
	Unit* temp = local::head;

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){
		
		// allocate
		int size_plane =(temp -> n + 1) * (temp -> m + 1);
		for(int i = 0; i < dg_fun::num_of_equation; ++i){
	
			result::error[i] = std::vector<double>(size_plane);
			result::exact[i] = std::vector<double>(size_plane);
		}

		// element size
		double del_x = (temp -> xcoords[1]) - (temp -> xcoords[0]); 
		double del_y = (temp -> ycoords[1]) - (temp -> ycoords[0]);  

		// wave ---------------------------------------------------------------------------------------		
		Exact_solution_Gaussian(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
					del_x, del_y, result::exact, dg_time::t_total);
//		Exact_solution_Gaussian2(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
//					del_x, del_y, result::exact, dg_time::t_total);
		// --------------------------------------------------------------------------------------------

		// test ---------------------------------------------------------------------------------------
//		Exact_solution_sin(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
//					del_x, del_y, result::exact, dg_time::t_total);
		//----------------------------------------------------------------------------------------------

		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

			int nodei{};
	
			for(int i = 0; i <= temp -> n; ++i){
		
				double weight_x = nodal::gl_weights[temp -> n][i];

				for(int j = 0; j <= temp -> m; ++j){

					result::error[equ][nodei] = std::abs(result::exact[equ][nodei] 
									- temp -> solution[equ][nodei]);
	
//					result::L2_norm[equ] += result::error[equ][nodei] * result::error[equ][nodei]
//								* nodal::gl_weights[temp -> m][j] * weight_x
//								* del_x * del_y / 4.0;
					result::L2_norm[equ] += result::error[equ][nodei] * result::error[equ][nodei];


					++nodei;
				}
			}

		}

		result::error.clear();
		result::exact.clear();

		temp = temp -> next;	

	}

	// sum up L2_norm on processor 0
	std::vector<double> L2_recv(dg_fun::num_of_equation);
	MPI_Reduce(&result::L2_norm[0], &L2_recv[0], dg_fun::num_of_equation, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(mpi::rank == 0){
		//std::cout.precision(17);
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
			result::L2_norm[equ] = sqrt(L2_recv[equ]);
//			std::cout<< std::fixed <<result::L2_norm[equ] << "\n";
			std::cout << result::L2_norm[equ] << "\n";
		}
	}

}

/// @brief
/// Serial write exact error to file. 
void Write_error(){

	for(int k = 0; k < mpi::num_proc; ++k){

		if(mpi::rank == k){
			Write_error_each_proc(dg_time::t_total);

		}
		else{

			MPI_Barrier(MPI_COMM_WORLD);		
		
		}
		

	}
	

}

/// @brief
/// Write error to files. L2 norm per element. 
void Write_error_each_proc(double t){

	// generate the file name
	std::string filename = 	fileinfo::exact_error_filename;
	std::ofstream myfile; 	// stream class to write on files	


	// processor open file
	if(mpi::rank == 0){
		
		myfile.open(filename, std::ios::trunc);	// truncate the old file

		// headers
		myfile<< "TITLE = \"EXACT ERROR\" \n";
		myfile<< "VARIABLES = \"X\", \"Y\", \"RANK\", \"HLEVEL\", \"PORDER\"," 
				"\"ERR_P\",\"ERR_U\", \"ERR_V\", \n";
		
	}
	else{

		myfile.open(filename, std::ios::app);	// All output operations are performed at the end of the file
	}


	std::unordered_map<int, std::vector<double>> exact;
	std::unordered_map<int, std::vector<double>> error;
	std::vector<double> L2_norm;

	// temperary pointer
	Unit* temp = local::head;

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){
		
		// allocate
		int size_plane =(temp -> n + 1) * (temp -> m + 1);
		for(int i = 0; i < dg_fun::num_of_equation; ++i){
	
			error[i] = std::vector<double>(size_plane);
			exact[i] = std::vector<double>(size_plane);
		}

		// element size
		double del_x = (temp -> xcoords[1]) - (temp -> xcoords[0]); 
		double del_y = (temp -> ycoords[1]) - (temp -> ycoords[0]);  

		// wave ---------------------------------------------------------------------------------------		
		Exact_solution_Gaussian(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
					del_x, del_y, exact, dg_time::t_total);
//		Exact_solution_Gaussian2(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
//					del_x, del_y, result::exact, dg_time::t_total);
		// --------------------------------------------------------------------------------------------

		// test ---------------------------------------------------------------------------------------
//		Exact_solution_sin(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
//					del_x, del_y, result::exact, dg_time::t_total);
		//----------------------------------------------------------------------------------------------

		// get the L2 norm of the current element ----------------------------------------------------
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
			L2_norm.push_back(0.0);

			int nodei{};

			for(int i = 0; i <= temp -> n; ++i){
		
				double weight_x = nodal::gl_weights[temp -> n][i];

				for(int j = 0; j <= temp -> m; ++j){
	
					error[equ][nodei] = std::abs(exact[equ][nodei] - 
								temp -> solution[equ][nodei]);

					L2_norm[equ] += error[equ][nodei] * error[equ][nodei] 
							* nodal::gl_weights[temp -> m][j] * weight_x
							* del_x * del_y / 4.0;

					++nodei;
				}
			}

			L2_norm[equ] = sqrt(L2_norm[equ]);

		}
		// --------------------------------------------------------------------------------------------


		// write result to file --------------------------------------------------
		myfile << std::fixed;
		myfile << std::setprecision(5);
		myfile << "ZONE T= " << "\"" << "IEL" << std::setw(8) << g_elem << "\"," << "  " 
			<<"I=2, J=2, "<< "SOLUTIONTIME=" << t <<", DATAPACKING = POINT" << "\n";
		
		++g_elem;
		

		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[0] 
			<< "  " << mpi::rank << "  " << temp -> index[2]
			<< "  "<< temp -> n << "  " << L2_norm[0] << "  " << L2_norm[1] << "  "<< L2_norm[2] <<"\n";


		myfile << temp -> xcoords[0] << "  " << temp -> ycoords[1] 
			<< "  " << mpi::rank << "  " << temp -> index[2]
			<< "  "<< temp -> n << "  " << L2_norm[0] << "  " << L2_norm[1] << "  "<< L2_norm[2] <<"\n";


		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[0] 
			<< "  " << mpi::rank << "  " << temp -> index[2]
			<< "  "<< temp -> n << "  " << L2_norm[0] << "  " << L2_norm[1] << "  "<< L2_norm[2] <<"\n";


		myfile << temp -> xcoords[1] << "  " << temp -> ycoords[1] 
			<< "  " << mpi::rank << "  " << temp -> index[2]
			<< "  "<< temp -> n << "  " << L2_norm[0] << "  " << L2_norm[1] << "  "<< L2_norm[2] <<"\n";
		
		//-----------------------------------------------------------------------


		error.clear();
		exact.clear();
		L2_norm.clear();

		temp = temp -> next;	

	}

	myfile.close();
}

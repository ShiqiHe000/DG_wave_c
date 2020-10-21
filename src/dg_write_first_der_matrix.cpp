#include <iostream>
#include "dg_write_first_der_matrix.h"
#include <mpi.h>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "dg_single_index.h"

/// @brief
/// Output the first derivative matrix into file. 
/// Only rank 0 does the job. 
/// Assuming polynomial order x == y. 
/// @param n polynomial order. 
void Write_first_der_matrix(int n){

	if(mpi::rank == 0){

		// generate the file name
		std::stringstream ss;
		ss << fileinfo::first_der_place <<"der" << std::setfill('0') << std::setw(2) << n << ".txt";
		std::string filename = 	ss.str();
		std::ofstream myfile; 	// stream class to write on files	

		// processor open file
		myfile.open(filename, std::ios::trunc);	// truncate the old file

		myfile<< std::fixed;
		myfile<< std::setprecision(17);

		int index{};

		for(int i = 0; i <= n; ++i){
			for(int j = 0; j <= n; ++j){

				index = Get_single_index(i, j, n + 1);

				if(j < n){

					myfile << nodal::first_der[n][index] << ",";
				}
				else{

					myfile << nodal::first_der[n][index];
				}
			}
			myfile << "\n";
		}

		myfile.close();
		

	}
}

void Write_mfirst_der_matrix(int n){

	if(mpi::rank == 0){

		// generate the file name
		std::stringstream ss;
		ss << fileinfo::first_der_place <<"mder" << std::setfill('0') << std::setw(2) << n << ".txt";
		std::string filename = 	ss.str();
		std::ofstream myfile; 	// stream class to write on files	

		// processor open file
		myfile.open(filename, std::ios::trunc);	// truncate the old file

		myfile<< std::fixed;
		myfile<< std::setprecision(5);

		int index{};

		for(int i = 0; i <= n; ++i){
			for(int j = 0; j <= n; ++j){

				index = Get_single_index(i, j, n + 1);

				if(j < n){

					myfile << nodal::mfirst_der[n][index] << ",";
				}
				else{

					myfile << nodal::mfirst_der[n][index];
				}
			}
			myfile << "\n";
		}

		myfile.close();
		

	}
}

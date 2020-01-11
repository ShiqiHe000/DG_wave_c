#include <mpi.h>
#include "dg_single_index.h"


/// @brief convert the 2d array index into a single 1d array index
/// @param row row number (start with 0) 
/// @param col column number (start with 0)
/// @param total_col total number of column (start with 1)
int Get_single_index(int row, int col, int total_col){
	
	return (row * total_col + col);

}

/// @brief convert the 3d array index (i, j, k) into a single 1d array index
/// @param i row number (start with 0) 
/// @param j column number (start with 0)
/// @param k thrid coordinate number (start with 0)
/// @param total_row total number of row (start with 1)
/// @param total_col total number of column (start with 1)
int Get_single_index_3d(int i, int j, int k, int total_row, int total_col){

	return(i * total_col + j + total_row * total_col * k);

}

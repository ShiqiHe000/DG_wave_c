#include <mpi.h>
#include "dg_single_index.h"


/// @brief convert the 2d array index into a single 1d array index
/// @param row row number (start with 0) 
/// @param col column number (start with 0)
/// @param total_col total number of column (start with 1)
int Get_single_index(int row, int col, int total_col){
	
	return (row * total_col + col);

}

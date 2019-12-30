#include <mpi.h>
#include <unordered_map>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_basis.h"

/// @brief
/// Processor's local storages.
/// Generate the basis information (GL PIOINTS, WEIGHTS, DERIVATIVE 
/// MATRIX, etc.) for all polynomial level at the beginning,
/// and store it in hash tables. 
/// Therefore, after refinement start we do not need to 
/// re-generate the basis information, or store duplicated data. 
void Construct_basis_storage(){

	// x direction
	

}


/// @brief
/// GL points and weights, first order derivative matrix
/// @param n polynomial order
/// @param gl_p GL points
/// @param gl_weight GL weights
/// @param first_der first order derivative matrix
void Get_nodal_2d_storage_basis(int n, double* gl_p, double* gl_weight, double* first_der){
	
	GL(n, gl_p, gl_weight);
	
}

#include <mpi.h>
#include "dg_nodal_2d_storage.h"
#include "dg_local_storage.h"
#include "dg_basis_storage.h"
#include "dg_constructor.h"

/// @brief 
/// A nodal DG 2D constructor.
/// This constructor computes the GL nodes and weights, Lagrange interpolating
/// polynomials at the boundaries, first order derivative matrices and 
/// modified them.
void Construct_basis(){
	
	// construct DG basis
	Construct_basis_storage();	
	
}

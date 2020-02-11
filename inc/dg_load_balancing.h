#ifndef DG_LOAD_BALANCING_H
#define DG_LOAD_BALANCING_H

#include <vector>

/// @brief
/// This structure assist to form the processor mapping table.
/// @param iproc Process's rank number.
/// @param gnum first element's global number.
struct pmap{

	int irank;	

	int gnum;

};

/// @brief
/// An envelope contains the element to be send to former proc and latter proc. 
/// @param pre elements whose destination is former proc. The vector stores the key of the atrget elements.  
/// @param next elements whose destination is the next proc. The vector stores the key of the atrget elements. 
struct sending_envelope{

	std::vector<int> pre;

	std::vector<int> next;

};

/// @brief
/// Data Structure to support to form the ownership table.
/// @param local_key key of local element who is on the MPI boundary.
/// @param owner_rank the destination rank this element will be after repatitioning. 
/// @param hlevel current element's h-refinement level. 
struct ownership{

	int local_key;
	
	int owner_rank;

	int hlevel;
};

// functions-----------------------------------------------------------------------------
void Build_mapping_table();

//void Form_ownership_table();

void Update_mpi_boundary();
// --------------------------------------------------------------------------------------

#endif

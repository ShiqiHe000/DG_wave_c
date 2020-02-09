#ifndef DG_LOAD_BALANCING_H
#define DG_LOAD_BALANCING_H

#include <vector>

/// @brief
/// This structure assist to form the processor mapping table.
/// @param iproc Process's rank number.
/// @param gnum first element's global number.
/// @param key first element's key.
struct pmap{

	int irank;	

	int gnum;

//	int key;

};

/// @brief
/// An envelope contains the element to be send to former proc and latter proc. 
/// @param pre elements whose destination is former proc. The vector stores the key of the atrget elements.  
/// @param next elements whose destination is the next proc. The vector stores the key of the atrget elements. 
struct sending_envelope{

	std::vector<int> pre;

	std::vector<int> next;

};


void Build_mapping_table();


#endif

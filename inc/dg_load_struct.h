#ifndef DG_LOAD_STRUCT_H
#define DG_LOAD_STRUCT_H

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

	std::vector<long long int> pre;

	std::vector<long long int> next;

};


#endif

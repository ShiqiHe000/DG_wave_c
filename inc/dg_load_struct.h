#ifndef DG_LOAD_STRUCT_H
#define DG_LOAD_STRUCT_H

#include <vector>

/// @brief
/// This structure assist to form the processor mapping table. 
/// The last variable prefix_sum is used to evaluate the optimal bottleneck. 
/// @param iproc Process's rank number.
/// @param gnum first element's global number.
/// @param prefix_sum the prefix sum of this element. 
//struct pmap_quality{
//
//	int irank;	
//
//	int gnum;
//
//	double prefix_sum;
//
//};

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

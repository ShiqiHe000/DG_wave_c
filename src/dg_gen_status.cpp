#include "dg_gen_status.h"
#include <cassert>
#include <cmath>

/// @brief 
/// Generate the Hilbert status.
/// Restriction: square domain
/// @param n pow(2, n) : element number on each boundary in exponential form.
/// @param status element status array.
void Gen_status(int n, char* status){
	
	// make sure element number is valid.
	assert(n >= 0 && "n must >= 0");
	
	// total element number
	const int total_elem = pow(2, 2 * n);

	const char four[]{'A', 'H', 'H', 'B', '\0'};
	
	// 3 situation, depending on the element number
	if(n > 1){	// more than 4 elem
		
		int div = total_elem / 4;
		int* find = new int[n-1]{};	// record path

		for(int k = 0; k < total_elem; ++k){
			int mid{};
			for(int s = 0; s < n-1; ++s){
				mid = k % div;
				find[s] = k / div;
				div /= 4;	
			}
			
			status[k] = four[find[n-2]];
			
			for(int s = 0; s < n-2 ; ++s){
				//status[k] = 

			}

		}
		

	}
	else if(n == 1){	// 4 elem
		for(int k = 0; k < total_elem; ++k){
			status[k] = four[k];

		}
		return;
	}
	else{	// 1 elem

		status[0] = 'H';
		return;
	}


	int div = total_elem / 4;
	while(div >= 4){
		

	}

}

#include "dg_status_table.h"
#include <unordered_map>
#include <array>
#include <cassert>


/// @brief
/// Hilbert status lookup table. Input current status and column number, output new status.
/// @param input current status.
/// @param i column number. 
char Status_table( char input, int i){
	
	// static hash table for status
	// key: current status
	// mapped_value: array of 4 element status. 	
	static std::unordered_map<char, std::array<char, 5>> status_lookup = 
				{{'H', {'A', 'H', 'H', 'B', '\0'}}, 
				 {'A', {'H', 'A', 'A', 'R', '\0'}}, 
				 {'R', {'B', 'R', 'R', 'A', '\0'}},
				 {'B', {'R', 'B', 'B', 'H', '\0'}}};
	
	// check i
	assert( i >= 0 && i <= 3 && "column number should between 0~3." );
	
	std::array<char, 5> m = status_lookup[input];
	char out = m[i];

	return out;	
	
}

#include "dg_unit.h"
#include "dg_param.h"
#include <vector>

//int Cantor_pairing(int x, int y){
//	
//	return ((x + y) * (x + y + 1) / 2 + y);
//}


/// @brief
/// constructor (default).
/// default: level = 0, n = min poly order, m = max poly order
/// index = {0, 0, 0}, mpi_f = false, x/ycoords = {0.0, 0.0}
/// solution = nullptr. 
Unit::Unit() : level(0), n(grid::nmin), m(grid::nmin)
{
	facen = std::vector<std::vector<Face>>(4);
}

/// @brief
/// Get the unique key based on the element index (i, j, k).
/// Return 
//int Unit::GetKey(){
//	
//	int key = Cantor_pairing(index[0], index[1]);
//	key = Cantor_pairing(index[2], key);
//
//	return key;
//
//}
//

/// @brief
/// Cantor pairing function: a bijection between n*n -> n. 
/// @param x number 1
/// @param y number 2
//int Unit::Cantor_pairing(int x, int y){
//	
//	return ((x + y) * (x + y + 1) / 2 + y);
//}




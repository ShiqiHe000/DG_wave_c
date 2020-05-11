#include "dg_cantor_pairing.h"
#include <cassert>

/// @brief 
/// Cantor pairing function. Bijection two integers into one unique integer.
/// @param x integer 1.
/// @param y integer 2.
long long int Cantor_pairing_fun(long long int x, long long int y){
	
	return ((x + y) * (x + y + 1) / 2 + y);

} 


/// @brief 
/// Pairing function, more spatial efficient than Cantor pairing.
///  Bijection two integers into one unique integer.
/// @param x integer 1.
/// @param y integer 2.
long long int Szudzik(long long int x, long long int y){

	assert(x >= 0 && y >= 0 && "The pairing function only deals with natural numbers. \n");

	long long int key = x >= y ? (x * x + x + y) : (y * y + x);

	return key;
}


/// @brief
/// Input element index, output the key
/// @param i x index. 
/// @param j y index. 
/// @param k level.  
long long int Get_key_fun(int i, int j, int k){
 
	// Cantor pairing function =============================================================	
//	long long int key = Cantor_pairing_fun((long long int)i, (long long int)j);
//	key = Cantor_pairing_fun(key, (long long int)k);
	// =====================================================================================


	// Szudzik pairing function ============================================================
	long long int key = Szudzik((long long int)i, (long long int)j);
	key = Szudzik(key, (long long int)k);
	// =====================================================================================

	return key;

}

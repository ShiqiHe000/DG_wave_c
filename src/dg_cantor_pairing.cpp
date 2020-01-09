#include "dg_cantor_pairing.h"
/// @brief 
/// Cantor pairing function. Bijection two integers into one unique integer.
/// @param x integer 1.
/// @param y integer 2.
int Cantor_pairing_fun(int x, int y){
	
	return ((x + y) * (x + y + 1) / 2 + y);

} 


/// @brief
/// Input element index, output the key
/// @param i x index. 
/// @param j y index. 
/// @param k level.  
int Get_key_fun(int i, int j, int k){
 	
	int key = Cantor_pairing_fun(i, j);
	key = Cantor_pairing_fun(key, k);

	return key;

}

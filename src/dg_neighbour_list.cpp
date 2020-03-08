#include "dg_neighbour_list.h"
#include <unordered_map>
#include <vector>
#include "dg_cantor_pairing.h"
#include "dg_param.h"
#include <cmath>

// forward declaration-----------------------------------------------------
void Total_neighbours(int level_now, int& all, int& my_position);

void Neighbours_array_x(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours);

void Neighbours_array_y(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours);
//-------------------------------------------------------------------------


/// @brief
/// Compute the total number of possible neighbours.
/// @param level_now The current element h-refinement level.
/// @param all Total position neighbour number.
/// @param my_position The same size neighbour's number. 
void Total_neighbours(int level_now, int& all, int& my_position){

	for(int i = grid::hlevel_max; i >= 0; --i){

		if(i > level_now){

			all += (int)std::pow(2, i - level_now);

		}
		else if(i == level_now){

			++all;
			my_position = all - 1;
		}
		else{
			++all;
		}
	}

} 

/// @brief
/// Form possible array in h-refinement level decending sequence (x direction).
/// @param i same size neighbour's x integer coordinate.
/// @param j same size neighbour's y integer coordinate.
/// @param k same size neighbour's h-refinement level.
/// @param local_key the key of the current element.
/// @param facen facen face number. 
/// @param neighbours neighbours list. 
void Neighbours_array_x(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours){

	int x{}; 
	
	if(facen == 0){	// south
		x = 1;
	}

	int total_n{};
	int same_size{};

	Total_neighbours(k, total_n, same_size);

	neighbours[local_key] = std::vector<int>(total_n);	// allocate space

	// same size neighbour
	neighbours[local_key][same_size] = Get_key_fun(i, j, k);

	std::vector<int> pre_level{i, j, k};

	int elem = same_size - 1;

	for(int m = k + 1; m <= grid::hlevel_max; ++m){

		int layer_elem = (m - k) * 2;
		
		pre_level[0] = 2 * pre_level[0] + x;
		pre_level[1] = 2 * pre_level[1];
		pre_level[2]++;

		for(int n = layer_elem - 1; n >= 0; --n){
			int key = Get_key_fun(pre_level[0], pre_level[1] + n, pre_level[2]);				
			
			neighbours[local_key][elem] = key;

			elem--;
		}
	}

	pre_level[0] = i / 2;
	pre_level[1] = j / 2;
	pre_level[2] = k - 1;

	elem = same_size + 1;

	for(int m = k - 1; m >= 0; --m){

		int key = Get_key_fun(pre_level[0], pre_level[1], pre_level[2]);

		neighbours[local_key][elem] = key;

		++elem;
	
		pre_level[0] /= 2;
		pre_level[1] /= 2;
		pre_level[2]--;
	}

}


/// @brief
/// Form possible array in h-refinement level decending sequence (y direction).
/// @param i same size neighbour's x integer coordinate.
/// @param j same size neighbour's y integer coordinate.
/// @param k same size neighbour's h-refinement level.
/// @param local_key the key of the current element.
/// @param facen facen face number. 
/// @param neighbours neighbours list. 
void Neighbours_array_y(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours){

	int y{};
	if(facen == 2){
		y = 1;
	}
	
	int total_n{};
	int same_size{};

	Total_neighbours(k, total_n, same_size);

	neighbours[local_key] = std::vector<int>(total_n);	// allocate space

	// same size neighbour
	neighbours[local_key][same_size] = Get_key_fun(i, j, k);

	std::vector<int> pre_level{i, j, k};

	int elem = same_size - 1;

	for(int m = k + 1; m <= grid::hlevel_max; ++m){

		int layer_elem = (m - k) * 2;
		
		pre_level[0] = 2 * pre_level[0];
		pre_level[1] = 2 * pre_level[1] + y;
		pre_level[2]++;

		for(int n = layer_elem - 1; n >= 0; --n){
			int key = Get_key_fun(pre_level[0] + n, pre_level[1], pre_level[2]);				
			
			neighbours[local_key][elem] = key;

			elem--;
		}
	}

	pre_level[0] = i / 2;
	pre_level[1] = j / 2;
	pre_level[2] = k - 1;

	elem = same_size + 1;

	for(int m = k - 1; m >= 0; --m){

		int key = Get_key_fun(pre_level[0], pre_level[1], pre_level[2]);

		neighbours[local_key][elem] = key;

		++elem;
	
		pre_level[0] /= 2;
		pre_level[1] /= 2;
		pre_level[2]--;
	}
}

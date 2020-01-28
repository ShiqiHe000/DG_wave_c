//#include "dg_adapt.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include <cstdlib>	// random number
#include "dg_cantor_pairing.h"	// for the key
#include <cassert>
#include <unordered_map>
#include "dg_status_table.h"
#include "dg_param.h"

/// global variable
double xcoord_new[3]{};
double ycoord_new[3]{};

/// forward declaration
void Get_coordinates(int ith, double* xcoord, double* ycoord);
void Gen_index(int ith, int* index, int* index_new);

/// @brief
/// Random h-refinement scheme. Each element has 30% chance to split.
void h_refinement(){

	Unit* temp = local::head;
	Unit* temp2 = temp;

	int increment{};
	
	for(int k = 0; k < local::local_elem_num; ++k){
		
		// generate random number
		int rand_num = rand() % 10 + 1;	// random number between [1, 10]

		bool check = ((temp -> index[2]) < grid::hlevel_max ) ? true : false;

		if(rand_num <= 3 && check){	// refine
			
			increment += 3; 
			
			// new coordinates
			xcoord_new[0] = temp -> xcoords[0]; xcoord_new[2] = temp -> xcoords[1];
			ycoord_new[0] = temp -> ycoords[0]; ycoord_new[2] = temp -> ycoords[1];
			xcoord_new[1] = ((temp -> xcoords[0]) + (temp -> xcoords[1])) / 2.0;
			ycoord_new[1] = ((temp -> ycoords[0]) + (temp -> ycoords[1])) / 2.0;

			// current key
			int old_key = Get_key_fun(temp->index[0], temp->index[1], temp->index[2]);
		
			// previous key between 4 children	
			int pre_key{};

			// create 4 unit in hash table
			for(int i = 0; i < 4; ++i){
				

				int position = Sibling_position(temp->status, i);

				int index_new[3]{};
				Gen_index(position, temp->index, index_new);
			
				int new_key = Get_key_fun(index_new[0], index_new[1], index_new[2]);

				// create unit
				local::Hash_elem[new_key] = new Unit();

				if(i == 0 && k == 0){	// first elem, update head
					local::head = local::Hash_elem[new_key];
					
				}
				else if(i == 0 && k != 0){
					temp2 -> next = local::Hash_elem[new_key];

				}
				
				// unit index
				local::Hash_elem[new_key] -> index[0] = index_new[0];
				local::Hash_elem[new_key] -> index[1] = index_new[1];
				local::Hash_elem[new_key] -> index[2] = index_new[2];

				// status
				local::Hash_elem[new_key] -> status = Status_table(temp->status, i );

				// coordinates
				Get_coordinates(position, local::Hash_elem[new_key] -> xcoords,
						 local::Hash_elem[new_key] -> ycoords);


				// ith child
				local::Hash_elem[new_key] -> child_position = position;

				// poly order
				local::Hash_elem[new_key] -> n = temp -> n;
				local::Hash_elem[new_key] -> m = temp -> m;
				
				// faces(here we construct between siblings)
				
					

				// interpolate solutions



				// form link
				if(i > 0){
					local::Hash_elem[pre_key] -> next = local::Hash_elem[new_key];
				}

				pre_key = new_key;
				
				temp2 = local::Hash_elem[new_key];
				
			}	
			
			temp2 -> next = temp -> next;

			// erase the parent
			local::Hash_elem.erase(old_key);
			temp = temp2;
		}
		temp2 = temp;
		temp = temp -> next;
	}

	local::local_elem_num += increment;

}

/// @brief
/// Build faces between siblngs. Each element has two faces adjecant to siblings. 
/// @param new_key current child element's key
/// @param position ith child. 
void Two_siblings(int new_key, int position){
	
	// four faces
	for(int i = 0; i < 4; ++i){

		bool t = Sibling_table(position, i);

		if(t){
			local::Hash_elem[new_key] -> facen[i].face_type = 'L';
			local::Hash_elem[new_key] -> facen[i].hlevel = local::Hash_elem[new_key] -> index[2];
			local::Hash_elem[new_key] -> facen[i].porderx = local::Hash_elem[new_key] -> n;
			local::Hash_elem[new_key] -> facen[i].pordery = local::Hash_elem[new_key] -> m;

			// key
			if(i == 0){	// south
				
				local::Hash_elem[new_key] -> facen[i].key = Get_key_fun(local::Hash_elem[new_key] -> index[0] - 1, 
											local::Hash_elem[new_key] -> index[1], 
											local::Hash_elem[new_key] -> index[2]);
		

			}
		}

	}

}


/// @brief
/// Form the element face vector
/// @param position Relative position of 4 children.
/// @param new_key child's key.
/// @param old_key parent's key.
void Get_facen(int position, int new_key, int old_key){

	assert(position >=0 && position <= 3 && "Error: The ith child's position is out of range.");


	if(position == 0){

//		for()
//		local::Hash_elem[new_key] -> facen[]		


	}
	else if(position == 1){


	}
	else if(position ==2){


	}
	else{



	}

}

/// @brief
/// Generate the new unit index based on the child position and its old unit index
/// @param ith ith child
/// @param index parent's unit index
/// @param index_new child's index (size 3)
void Gen_index(int ith, int* index, int* index_new){

	assert(ith >=0 && ith<= 3 && "Error: The ith child's position is out of range.");
	
	if(ith == 0){

		index_new[0] = 2 * index[0];
		index_new[1] = 2 * index[1];
		index_new[2] = index[2] + 1;
	}
	else if(ith == 1){

		index_new[0] = 2 * index[0] + 1;
		index_new[1] = 2 * index[1];
		index_new[2] = index[2] + 1;

	}
	else if(ith == 2){

		index_new[0] = 2 * index[0] + 1;
		index_new[1] = 2 * index[1] + 1;
		index_new[2] = index[2] + 1;

	}
	else{


		index_new[0] = 2 * index[0];
		index_new[1] = 2 * index[1] + 1;
		index_new[2] = index[2] + 1;

	}

}

/// @brief 
/// Generate the new coordinates base on the child position inside the quard
/// @param ith ith child
/// @param xcoord new unit's xcoord
/// @param ycoord new unit's ycoord
void Get_coordinates(int ith, double* xcoord, double* ycoord){

	assert(ith >=0 && ith<= 3 && "Error: The ith child's position is out of range.");

	if(ith == 0){
		xcoord[0] = xcoord_new[0]; xcoord[1] = xcoord_new[1];
		ycoord[0] = ycoord_new[0]; ycoord[1] = ycoord_new[1];
	}
	else if(ith == 1){

		xcoord[0] = xcoord_new[1]; xcoord[1] = xcoord_new[2];
		ycoord[0] = ycoord_new[0]; ycoord[1] = ycoord_new[1];

	}
	else if(ith == 2){

		xcoord[0] = xcoord_new[1]; xcoord[1] = xcoord_new[2];
		ycoord[0] = ycoord_new[1]; ycoord[1] = ycoord_new[2];
	}
	else{

		xcoord[0] = xcoord_new[0]; xcoord[1] = xcoord_new[1];
		ycoord[0] = ycoord_new[1]; ycoord[1] = ycoord_new[2];
	}

}

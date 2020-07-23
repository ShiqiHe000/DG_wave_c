#include "dg_amr.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include <cstdlib>	// random number
#include "dg_cantor_pairing.h"	// for the key
#include <cassert>
#include <unordered_map>
#include "dg_status_table.h"
#include "dg_param.h"
#include <array>
#include <cassert>
#include "dg_elem_length.h"
#include "dg_neighbour_list.h"
#include <algorithm>
#include "dg_interpolate_to_new_points.h"
#include "dg_prefinement.h"
#include <iostream>	// test
#include "dg_test.h"	// test
#include "dg_single_index.h"	// test

/// global variable
double xcoord_new[3]{};
double ycoord_new[3]{};

/// forward declaration ----------------------------------------------------------------------------
void hpc_refinement();

void Get_coordinates(int ith, double* xcoord, double* ycoord);

void Gen_index(int ith, int* index, int* index_new);

void Two_siblings(long long int new_key, int position);

void Four_children(int ith, int i, int j, int level, std::array<long long int, 4>& children);

void Non_sibling_interfaces(Unit* last, long long int parent, std::array<long long int, 4>& children);

void Form_one_direction(long long int key1, long long int key2, long long int parent, int facen);

void Match_neighbours(long long int parent, long long int local_key, int facen, std::vector<long long int>& neighbours);

void Flag_elem(int kt);

void Coarsen_criteria(Unit* temp, bool& pass, std::array<long long int, 4>& four_keys);

int Parent_position(int i, int j);

void Inherit_from_children(long long int c1, long long int c2, long long int p_key, int facen);

void Form_parent_faces(std::array<long long int, 4>& four_keys, long long int p_key);

void Change_neighbour_coasen_case1(long long int c1, int facen, long long int p_key);

void Change_neighbour_coarsen_case2(long long int c1, long long int c2, int facen, long long int p_key);

bool First_child(Unit* temp);

void Ref_coods_x(long long int key, int position, Unit* temp);
void Ref_coods_y(long long int key, int position, Unit* temp);

// void Coarsen_results(Unit* temp, bool pass);	// test
void Print_inter(long long int new_key, long long int old_key);	// test
// --------------------------------------------------------------------------------------------------

/// @brief
/// loop through all the elements and flag then for refining or coarsening. 
void Flag_elem(int kt){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		// generate random number
//		int rand_num = rand() % 10 + 1;	// random number between [1, 10]

//int key_test = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
//if(mpi::rank == 0 && kt == 0 ){
//
//	rand_num = 1;
//	
////	temp -> prefine = true;
//
//}
//else if(mpi::rank == 0 && kt != 0){
//
//	rand_num = 4;
//	
////	temp -> coarsen = true;
//}
//else{
//	rand_num = 4;
//}
int rand_num = 4;

		bool check_h = ((temp -> index[2]) < grid::hlevel_max ) ? true : false;
		bool check_c = ((temp -> index[2]) > 0)	 ? true : false;

		if(check_h && rand_num <= 3 ){	// h-refinement

			temp -> hrefine = true;
		}
		else if(check_c && rand_num > 5 ){	// coarse
			
			temp -> coarsen = true;
		}
		temp = temp -> next;
	}
}



/// @brief
/// Random h-refinement scheme. Each element has 30% chance to split.
void hpc_refinement(){

	Unit* temp = local::head;
	Unit* temp2 = temp;

	int increment{};	// new added element because of hrefinemnt
	int decrement{};	// element decrement because coarening
	
	for(int k = 0; k < local::local_elem_num; ++k){

		if(temp -> hrefine){	// h - refine======================================================================
			
			increment += 3; 
			
			// new coordinates
			xcoord_new[0] = temp -> xcoords[0]; xcoord_new[2] = temp -> xcoords[1];
			ycoord_new[0] = temp -> ycoords[0]; ycoord_new[2] = temp -> ycoords[1];
			xcoord_new[1] = ((temp -> xcoords[0]) + (temp -> xcoords[1])) / 2.0;
			ycoord_new[1] = ((temp -> ycoords[0]) + (temp -> ycoords[1])) / 2.0;
			
			// current key
			long long int old_key = Get_key_fun(temp->index[0], temp->index[1], temp->index[2]);
		
			// previous key between 4 children	
			long long int pre_key{};

			// create 4 unit in hash table
			for(int i = 0; i < 4; ++i){
				

				int position = Sibling_position(temp->status, i);
				
				int index_new[3]{};
				Gen_index(position, temp->index, index_new);
				
				long long int new_key = Get_key_fun(index_new[0], index_new[1], index_new[2]);

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
				
				// reference coods
				Ref_coods_x(new_key, position, temp);
				Ref_coods_y(new_key, position, temp);

				// poly order
				local::Hash_elem[new_key] -> n = temp -> n;
				local::Hash_elem[new_key] -> m = temp -> m;
				
				// faces(here we construct between siblings)
				Two_siblings(new_key, position);
					
				// form link
				if(i > 0){
					local::Hash_elem[pre_key] -> next = local::Hash_elem[new_key];
				}

				pre_key = new_key;
				
				temp2 = local::Hash_elem[new_key];
				
			}	
			// form then external interfaces between 4 children and updates their neighbour's faces			
			std::array<long long int, 4> children;	// four siblings' keys
			Non_sibling_interfaces(temp2, old_key, children);
			// interpolate solution to four children
			Solutions_to_children(children, old_key);

			temp2 -> next = temp -> next;

			// erase the parent
			local::Hash_elem.erase(old_key);
			temp = temp2;
		}
		else if(temp -> prefine){	// p-refinement===============================================================

			p_refinement_apply(temp);

		}	
		else if(temp -> coarsen){	// coarsening ===============================================================
			
			if(temp -> n > grid::nmin){	// if polynomial order is larger then the minimum
		
				// p-coarsening
//				p_coarsening_interpolate(temp);	// interpolation

				p_coarsening_L2(temp);	// L2 projection
			
			}
			else{	// h-coasrening

				if(temp -> index[2] > 0){	// only non-root element can be coarsen. 

					// only flaged element at the child posiiton 0 or 2 can croasen
					bool first = First_child(temp);
	
					if(first){
						std::array<long long int, 4> four_keys;
						bool pass; 	
						Coarsen_criteria(temp, pass, four_keys);
	
	//					Coarsen_results(temp, pass);
	
						if(pass){ // if four siblings all want to coarse
							decrement += 3;	
		
							// generate parent's key
							long long int key_p = Get_key_fun((temp->index[0]) / 2, 
									(temp->index[1]) / 2, temp->index[2] - 1);
							
							local::Hash_elem[key_p] = new Unit();	// create parent
						
							// index	
							local::Hash_elem[key_p] -> index[0] = (temp -> index[0]) / 2;
							local::Hash_elem[key_p] -> index[1] = (temp -> index[1]) / 2;
							local::Hash_elem[key_p] -> index[2] = (temp -> index[2]) - 1;
		
							// status
							local::Hash_elem[key_p] -> status = Go_back_to_parent(temp -> status);
	
							// coordinates (inherit from elements at position 0 and 2)
							local::Hash_elem[key_p]	-> xcoords[0] = 
										local::Hash_elem[four_keys[0]] -> xcoords[0];
	
							local::Hash_elem[key_p]	-> ycoords[0] = 
										local::Hash_elem[four_keys[0]] -> ycoords[0];
	
							local::Hash_elem[key_p]	-> xcoords[1] = 
										local::Hash_elem[four_keys[2]] -> xcoords[1];
	
							local::Hash_elem[key_p]	-> ycoords[1] = 
										local::Hash_elem[four_keys[2]] -> ycoords[1];
	
							// relative position
							local::Hash_elem[key_p] -> child_position = 
									Parent_position(local::Hash_elem[key_p] -> index[0],
											local::Hash_elem[key_p] -> index[1]);
							// poly orders
							local::Hash_elem[key_p] -> n = temp -> n; // the sibling should share same n, m
							local::Hash_elem[key_p] -> m = temp -> m; 
	
	
							// element boundary in reference space
							local::Hash_elem[key_p] -> ref_x[0] = 
										local::Hash_elem[four_keys[0]] -> ref_x[0];

							local::Hash_elem[key_p] -> ref_x[1] = 
										local::Hash_elem[four_keys[1]] -> ref_x[1];
	
							local::Hash_elem[key_p] -> ref_y[0] = 
										local::Hash_elem[four_keys[0]] -> ref_y[0];

							local::Hash_elem[key_p] -> ref_y[1] = 
										local::Hash_elem[four_keys[3]] -> ref_y[1];
	
							// adjust linked list
							if(k == 0){	// first elem
								local::head = local::Hash_elem[key_p];	
							}
							else{
	
								temp2 -> next = local::Hash_elem[key_p];
							}
							Unit* temp3 = temp;
							for(int m = 0; m < 3; ++m){	// rearch last sibling
	
								temp3 = temp3 -> next;
							}
							local::Hash_elem[key_p] -> next = temp3 -> next;
							
							temp = local::Hash_elem[key_p];	// move pointer to the last 
						
							k += 3;	// skip other siblings 
	
							// L2 project back to parent 						
							Solution_back_to_parent(four_keys, key_p);
	
							// form the face info + change neighbours faces
							Form_parent_faces(four_keys, key_p);
	
							// erase four siblings
							for(int i = 0; i < 4; ++i){
								local::Hash_elem.erase(four_keys[i]);
							}
						}
					}
				}
			}
			
		}
		temp2 = temp;
		temp = temp -> next;

	}

	local::local_elem_num += increment;
	local::local_elem_num -= decrement;

	
}

//// test function, trivial
//void Coarsen_results(Unit* temp, bool pass){
//
//	if(pass){
//
//		std::cout<< "rank "<< mpi::rank<<" coord " << temp -> index[0] << temp -> index[1]<< temp -> index[2]<< 
//			" state "<< temp -> status << " child_position "<< temp -> child_position <<"\n";
//	}
//
//}

/// @brief
/// Generate the reference boundary coordinate in x direciton.
/// @param key current child's key.
/// @param temp pointer to the parent element.
void Ref_coods_x(long long int key, int position, Unit* temp){

	assert(position >= 0 && position <= 3 && "position should be inside [0, 3]");

	double avg = ((temp -> ref_x[0]) + (temp -> ref_x[1])) / 2.0;

	if(position == 1 || position == 2){
		
		local::Hash_elem[key] -> ref_x[0] = avg;
		local::Hash_elem[key] -> ref_x[1] = temp -> ref_x[1];
	
	}
	else{	//position == 0 || position == 3

		local::Hash_elem[key] -> ref_x[0] = temp -> ref_x[0];
		local::Hash_elem[key] -> ref_x[1] = avg;

	}
}

/// @brief
/// Generate the reference boundary coordinate in y direciton.
/// @param key current child's key.
/// @param temp pointer to the parent element.
void Ref_coods_y(long long int key, int position, Unit* temp){

	assert(position >= 0 && position <= 3 && "position should be inside [0, 3]");

	double avg = ((temp -> ref_y[0]) + (temp -> ref_y[1])) / 2.0;

	if(position == 2 || position == 3){
		
		local::Hash_elem[key] -> ref_y[0] = avg;
		local::Hash_elem[key] -> ref_y[1] = temp -> ref_y[1];
	
	}
	else{	//position == 0 || position == 1

		local::Hash_elem[key] -> ref_y[0] = temp -> ref_y[0];
		local::Hash_elem[key] -> ref_y[1] = avg;

	}
}

/// @brief 
/// Record the neighbours info of new parent (4 sides). 
/// @param four_keys An array which stores the keys of four siblings in the sequence of relative posiiton. 
/// @param p_key the key of parent. 
void Form_parent_faces(std::array<long long int, 4>& four_keys, long long int p_key){

	// south
	Inherit_from_children(four_keys[0], four_keys[3], p_key, 0);

	// north
	Inherit_from_children(four_keys[1], four_keys[2], p_key, 1);

	// west
	Inherit_from_children(four_keys[0], four_keys[1], p_key, 2);

	// east
	Inherit_from_children(four_keys[3], four_keys[2], p_key, 3);
	
}


/// @brief
/// Examine whether this child is the first sibling.
/// @param temp pointer to the current child.
bool First_child(Unit* temp){

	if(temp -> child_position == 0){

		if(temp -> status == 'A' || temp -> status == 'H'){
			return true;
		}
		else{ return false;}
	}
	else if(temp -> child_position == 2){

		if(temp -> status == 'B' || temp -> status == 'R'){

			return true;
		}
		else{ return false;}
	}
	else{ return false;}
}


/// @brief
/// Form the face info for parent element in one direction (N, S, E, W).
/// @param c1 key of the 1st element of the corresponding direction. 
/// @param c2 key of the 2st element of the corresponding direction. 
/// @param p_key parent's key. 
/// @param facen face direction (0, 1, 2, 3).
void Inherit_from_children(long long int c1, long long int c2, long long int p_key, int facen){

	if(local::Hash_elem[c1] -> facen[facen].front().face_type == 'B'){	// if on the physical boundary

		local::Hash_elem[p_key] -> facen[facen].emplace_back();
		local::Hash_elem[p_key] -> facen[facen].back().face_type = 'B';

		return;
	}

	// if not physical boundary. Compare children and their neighbour's size. 
	// If same size, record both. If larger, then they are facing the same one, record one. 
	int c_level = local::Hash_elem[c1] -> index[2];
	int n_level = local::Hash_elem[c1] -> facen[facen].front().hlevel;

	if(c_level > n_level){	// neighbour is larger than children. erase two children + insert parent

		// copy the info
		local::Hash_elem[p_key] -> facen[facen] = local::Hash_elem[c1] -> facen[facen];
		
		if(local::Hash_elem[c1] -> facen[facen].front().face_type == 'L'){	// change neighbours

			Change_neighbour_coarsen_case2(c1, c2, facen, p_key);

		}

	}
	else{	// neighbours are smaller or equaled in size (compare with children). Erase children + insert parent.

		local::Hash_elem[p_key]	-> facen[facen] = local::Hash_elem[c1] -> facen[facen];

		local::Hash_elem[p_key] -> facen[facen].insert(local::Hash_elem[p_key] -> facen[facen].end(), 
								local::Hash_elem[c2] -> facen[facen].begin(),
								local::Hash_elem[c2] -> facen[facen].end());

		
		Change_neighbour_coasen_case1(c1, facen, p_key);
		Change_neighbour_coasen_case1(c2, facen, p_key);

	}


}

/// @brief
/// Change neighbours inferfaces info. Case one: neighbours size equals or smaller than children's size. 
/// @param c1 Key of child.
/// @param facen Children's facen number. 
/// @param p_key Parent's key.
void Change_neighbour_coasen_case1(long long int c1, int facen, long long int p_key){

	int n_dir = Opposite_dir(facen);	// neighbour face direction

	for(auto it = local::Hash_elem[c1] -> facen[facen].begin();
		it != local::Hash_elem[c1] -> facen[facen].end(); ++it){

		if(it -> face_type == 'L'){

			long long int n_key = it -> key;

			--local::Hash_elem[n_key] -> facen[n_dir].front().hlevel;	// hlevel - 1
		
			local::Hash_elem[n_key] -> facen[n_dir].front().porderx = local::Hash_elem[p_key] -> n;
			local::Hash_elem[n_key] -> facen[n_dir].front().pordery = local::Hash_elem[p_key] -> m;
			
			local::Hash_elem[n_key] -> facen[n_dir].front().key = p_key;

			local::Hash_elem[n_key] -> facen[n_dir].front().ref_x = local::Hash_elem[p_key] -> ref_x;
			local::Hash_elem[n_key] -> facen[n_dir].front().ref_y = local::Hash_elem[p_key] -> ref_y;
			
		}

	}



}


/// @brief
/// Change neighbours inferfaces info. Case two: neighbours size is larger than children's size. 
/// @param c1 Key of child 1.
/// @param c2 Key of child 2.
/// @param facen Children's facen number. 
/// @param p_key Parent's key.
void Change_neighbour_coarsen_case2(long long int c1, long long int c2, int facen, long long int p_key){

	int n_dir = Opposite_dir(facen);	// neighbour face direction

	long long int n_key = local::Hash_elem[c1] -> facen[facen].front().key;

	for(auto it = local::Hash_elem[n_key] -> facen[n_dir].begin(); 
		it != local::Hash_elem[n_key] -> facen[n_dir].end(); ){


		if(it -> key == c1 || it -> key == c2){

			it = local::Hash_elem[n_key] -> facen[n_dir].erase(it);
		}
		else{
			++it;
		}
	}

	// append parent's info at the end
	local::Hash_elem[n_key] -> facen[n_dir].emplace_back('L', local::Hash_elem[p_key] -> index[2], 
								local::Hash_elem[p_key] -> n, local::Hash_elem[p_key] -> m, 
								p_key, mpi::rank, local::Hash_elem[p_key] -> ref_x, 
								local::Hash_elem[p_key] -> ref_y);
}
/// @brief
/// Generate parent relative position base on the position of the first child.
/// @brief i integer coordinate in x direction.
/// @brief j integer coordinate in y direction.
int Parent_position(int i, int j){

	int a = i % 2;
	int b = j % 2;

	if(a == 0 && b == 0){
		return 0;

	}
	else if(a == 0 && b != 0){

		return 3;
	}
	else if(a != 0 && b == 0){

		return 1;
	}
	else{	// a, b != 0

		return 2;
	}
}


/// @brief
/// Enable coarsening only if 4 siblings all agree to coarse, i.e., their are all flagged to be coaren, and they share
/// same polynomial order and h-refinement level. 
/// @param temp unit pointer to the current element.
/// @param pass boolean variable. If true then coarse the grid. 
/// @param four_keys array that stores the four siblings key (in relative position order 0 1 2 3).
void Coarsen_criteria(Unit* temp, bool& pass, std::array<long long int, 4>& four_keys){

	assert(temp -> child_position == 0 || temp -> child_position == 2 && "children position prob");

	// only two circumstances
	if(temp -> child_position == 0){	
		
		int my_level = temp -> index[2];
		int my_order = temp -> n;

		// first check north and east
		if(temp -> facen[1].front().face_type != 'L') { // only neighbour are resided locally
			pass = false;
			return;
		}
		if(my_level != (temp -> facen[1].front().hlevel)) { // only my neighbour is the same size
			pass = false;
			return;
		}	
		if(my_order != (temp -> facen[1].front().porderx)){	// only neighbour is the same order
									// right now assume element porderx == pordery

			pass = false;
			return;
		}

		if(temp -> facen[3].front().face_type != 'L'){
			pass = false;
			return;
		}
		if(my_level != (temp -> facen[3].front().hlevel)){
			pass = false;
			return;
		}
		if(my_order != (temp -> facen[3].front().porderx)){	// only neighbour is the same order
									// right now assume element porderx == pordery

			pass = false;
			return;
		}
		
		// check the last sibling
		long long int north_key = temp -> facen[1].front().key;
		if(local::Hash_elem[north_key] -> facen[3].front().face_type != 'L'){
			pass = false;
			return;
		}
		if(local::Hash_elem[north_key] -> facen[3].front().hlevel != my_level){
			pass = false;
			return;
		}
		if(local::Hash_elem[north_key] -> facen[3].front().porderx != my_order){
			pass = false;
			return;
		}
		
		// generate four childrens key (0 1 2 3)
		four_keys[0] = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
		four_keys[1] = Get_key_fun(temp -> index[0] + 1, temp -> index[1], temp -> index[2]);
		four_keys[2] = Get_key_fun(temp -> index[0] + 1, temp -> index[1] + 1, temp -> index[2]);
		four_keys[3] = Get_key_fun(temp -> index[0], temp -> index[1] + 1, temp -> index[2]);

	}
	else{	// child_position == 2

		int my_level = temp -> index[2];
		int my_order = temp -> n;

		// first check south and west
		if(temp -> facen[0].front().face_type != 'L'){
			pass = false;
			return;
		}
		if(my_level != (temp -> facen[0].front().hlevel)){
			pass = false;
			return;
		}
		if(my_order != (temp -> facen[0].front().porderx)){
			pass = false;
			return;
		}

		if(temp -> facen[2].front().face_type != 'L'){
			pass = false;
			return;
		}
		if(my_level != (temp -> facen[2].front().hlevel)){
			pass = false;
			return;
		}
		if(my_order != (temp -> facen[2].front().porderx)){
			pass = false;
			return;
		}
		
		// check the last sibling
		long long int south_key = temp -> facen[0].front().key;
		if(local::Hash_elem[south_key] -> facen[2].front().face_type != 'L'){
			pass = false;
			return;
		}
		if(local::Hash_elem[south_key] -> facen[2].front().hlevel != my_level){
			pass = false;
			return;
		}
		if(local::Hash_elem[south_key] -> facen[2].front().porderx != my_order){
			pass = false;
			return;
		}

		// generate four childrens key (0 1 2 3)
		four_keys[0] = Get_key_fun(temp -> index[0] - 1, temp -> index[1] - 1, temp -> index[2]);
		four_keys[1] = Get_key_fun(temp -> index[0], temp -> index[1] - 1, temp -> index[2]);
		four_keys[2] = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
		four_keys[3] = Get_key_fun(temp -> index[0] - 1, temp -> index[1], temp -> index[2]);
	}

	// pass test
	pass = true;
}




/// @brief
/// Form the non-sibling interfaces of 4 children. 
/// @param last pointer to the 3th child.
/// @param parent parent key.
/// @param children four children's keys. 
void Non_sibling_interfaces(Unit* last, long long int parent, std::array<long long int, 4>& children){

	// get four children's keys
	Four_children(last -> child_position, last -> index[0], last -> index[1], last -> index[2], children);

	// south
	Form_one_direction(children[0], children[3], parent, 0);
	
	// north
	Form_one_direction(children[1], children[2], parent, 1);
		
	// west
	Form_one_direction(children[0], children[1], parent, 2);
	
	// east
	Form_one_direction(children[3], children[2], parent, 3);
}

/// @brief
/// Loop through parent's face and select out the nieghbours to the two children on one face.
/// Plus update the neighbours face info.
/// @param parent parent's key.
/// @param local_key child's key.
/// @parem facen face number.
/// @param neighbours vector of the possible neighbours.
void Match_neighbours(long long int parent, long long int local_key, int facen, std::vector<long long int>& neighbours){

	int l_tol{};
	int local_length = Elem_length(local::Hash_elem[local_key] -> index[2]);

	for(auto it = local::Hash_elem[parent] -> facen[facen].begin(); 
		it != local::Hash_elem[parent] -> facen[facen].end(); ++it){

		auto it_n = std::find(neighbours.begin(), neighbours.end(), it -> key);

		if(it_n != neighbours.end()){	// if find, stores
			
			local::Hash_elem[local_key] -> facen[facen].emplace_back(it);

			// update neighbours face info
			if(it -> face_type == 'L'){	// only local element
				
				int n_dir = Opposite_dir(facen);	// neighbour face direction
				long long int n_key = it -> key;	// neighbour's key

				// erase parent info, if exist
				for(auto it_nf = local::Hash_elem[n_key] -> facen[n_dir].begin();
					it_nf != local::Hash_elem[n_key] -> facen[n_dir].end(); ++it_nf){

					if(it_nf -> key == parent){	// find old parent info
				
						it_nf = local::Hash_elem[n_key] -> facen[n_dir].erase(it_nf);
						break;
					}

				}
				// child info
				local::Hash_elem[n_key] -> facen[n_dir].emplace_back('L', 
									local::Hash_elem[local_key] -> index[2], 
									local::Hash_elem[local_key] -> n, 
									local::Hash_elem[local_key] -> m, 
									local_key, mpi::rank,
									local::Hash_elem[local_key] -> ref_x,
									local::Hash_elem[local_key] -> ref_y);

			}
			

			l_tol += Elem_length(it -> hlevel);

			while(l_tol >= local_length){ break; }
		}
	
	}
	
}

/// @brief 
/// Form one direction's non-siblnig interfaces of two children.
/// @param key1 1st child's key
/// @param key2 2nd child's key (two children in index i/j ascending mode).
/// @param parent parent's key.
/// @param facen face direction. 
void Form_one_direction(long long int key1, long long int key2, long long int parent, int facen){
	
	// if on the physical boundary, only inherit. 
	if(local::Hash_elem[parent] -> facen[facen].front().face_type == 'B'){

		local::Hash_elem[key1] -> facen[facen].push_back(Unit::Face());
		local::Hash_elem[key1] -> facen[facen].front().face_type = 'B';
	
		local::Hash_elem[key2] -> facen[facen].push_back(Unit::Face());
		local::Hash_elem[key2] -> facen[facen].front().face_type = 'B';

		return; 
	}

	std::array<long long int, 2> two = {key1, key2};
	
	std::unordered_map<long long int, std::vector<long long int>> neighbours;

	// not on the physical boundary, form possible nieghbour list and search
	if(facen == 0){	// south

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_x(i - 1, j, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}
	}
	else if(facen == 1){	// north 

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_x(i + 1, j, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}

	}
	else if(facen == 2){	// west

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_y(i, j - 1, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}

	}
	else{	// east

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_y(i, j + 1, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}

	}

}


/// @brief
/// Input the last child's element index, then returns the array 
/// of 4 children's key in ascending child_position sequence. 
/// The children key sequence is: SW, NW, NE, SE (irrelevent of curve orientation.)
/// @param ith ith children (only support 1 and 3).
/// @param i integer coordinate in x direction.
/// @param j integer coordinate in y direction.
/// @param level element hlevel.
/// @param children children key array. 
void Four_children(int ith, int i, int j, int level, std::array<long long int, 4>& children){

	assert((ith == 1 || ith == 3) && "Child_position must be 1 or 3");

	if(ith == 1){	// 1st child
	
		children[0] = Get_key_fun(i - 1, j, level);
		children[1] = Get_key_fun(i, j, level);
		children[2] = Get_key_fun(i, j + 1, level);
		children[3] = Get_key_fun(i - 1, j + 1, level);

	}
	else{
		children[0] = Get_key_fun(i , j - 1, level);
		children[1] = Get_key_fun(i + 1, j - 1, level);
		children[2] = Get_key_fun(i + 1, j, level);
		children[3] = Get_key_fun(i, j, level);
		
	}


}


/// @brief
/// Build faces between siblngs. Each element has two faces adjecant to siblings. 
/// @param new_key current child element's key
/// @param position ith child. 
void Two_siblings(long long int new_key, int position){
	
	// four faces
	for(int i = 0; i < 4; ++i){

		bool t = Sibling_table(position, i);

		if(t){
			local::Hash_elem[new_key] -> facen[i].push_back(Unit::Face());
			local::Hash_elem[new_key] -> facen[i][0].face_type = 'L';
			local::Hash_elem[new_key] -> facen[i][0].hlevel = local::Hash_elem[new_key] -> index[2];
			local::Hash_elem[new_key] -> facen[i][0].porderx = local::Hash_elem[new_key] -> n;
			local::Hash_elem[new_key] -> facen[i][0].pordery = local::Hash_elem[new_key] -> m;
			local::Hash_elem[new_key] -> facen[i][0].rank = mpi::rank;
			
			// key
			if(i == 0){	// south
				
				local::Hash_elem[new_key] -> facen[i][0].key = Get_key_fun(
										local::Hash_elem[new_key] -> index[0] - 1, 
											local::Hash_elem[new_key] -> index[1], 
											local::Hash_elem[new_key] -> index[2]);

				double x0 = local::Hash_elem[new_key] -> ref_x[0];
				double x1 = local::Hash_elem[new_key] -> ref_x[1];
				double del_x = x1 -x0;
				
				local::Hash_elem[new_key] -> facen[i][0].ref_x[0] = x0 - del_x;
				local::Hash_elem[new_key] -> facen[i][0].ref_x[1] = x0;
				local::Hash_elem[new_key] -> facen[i][0].ref_y = local::Hash_elem[new_key] -> ref_y;
		

			}
			else if(i == 1){	// north
				
				local::Hash_elem[new_key] -> facen[i][0].key = 
									Get_key_fun(local::Hash_elem[new_key] -> index[0] + 1, 
											local::Hash_elem[new_key] -> index[1], 
											local::Hash_elem[new_key] -> index[2]);
		
				double x0 = local::Hash_elem[new_key] -> ref_x[0];
				double x1 = local::Hash_elem[new_key] -> ref_x[1];
				double del_x = x1 - x0;
				
				local::Hash_elem[new_key] -> facen[i][0].ref_x[0] = x1;
				local::Hash_elem[new_key] -> facen[i][0].ref_x[1] = x1 + del_x;
				local::Hash_elem[new_key] -> facen[i][0].ref_y = local::Hash_elem[new_key] -> ref_y;

			}
			else if(i == 2){	// west
				
				local::Hash_elem[new_key] -> facen[i][0].key = Get_key_fun(local::Hash_elem[new_key] -> index[0],
										local::Hash_elem[new_key] -> index[1] - 1, 
										local::Hash_elem[new_key] -> index[2]);
		
				double y0 = local::Hash_elem[new_key] -> ref_y[0];
				double y1 = local::Hash_elem[new_key] -> ref_y[1];
				double del_y = y1 - y0;
				
				local::Hash_elem[new_key] -> facen[i][0].ref_x = local::Hash_elem[new_key] -> ref_x;
				local::Hash_elem[new_key] -> facen[i][0].ref_y[0] = y0 - del_y;
				local::Hash_elem[new_key] -> facen[i][0].ref_y[1] = y0;

			}
			else{	// east
				
				local::Hash_elem[new_key] -> facen[i][0].key = Get_key_fun(local::Hash_elem[new_key] -> index[0],
										local::Hash_elem[new_key] -> index[1] + 1, 
										local::Hash_elem[new_key] -> index[2]);
		
				double y0 = local::Hash_elem[new_key] -> ref_y[0];
				double y1 = local::Hash_elem[new_key] -> ref_y[1];
				double del_y = y1 - y0;
				
				local::Hash_elem[new_key] -> facen[i][0].ref_x = local::Hash_elem[new_key] -> ref_x;
				local::Hash_elem[new_key] -> facen[i][0].ref_y[0] = y1;
				local::Hash_elem[new_key] -> facen[i][0].ref_y[1] = y1 + del_y;

			}
		}

	}

}



/// @brief
/// Generate the new unit index based on the child position and its old (parent) unit index
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

/// @brief
/// test function
/// Print out the interpolants
/// @parma new_key child's key
/// @param old_key parent's key
void Print_inter(long long int new_key, long long int old_key){

	Unit* temp_p = local::Hash_elem[old_key];
	Unit* temp_c = local::Hash_elem[new_key];

	// parent's solution
	std::cout<< "===========parent " << old_key << " child " << new_key <<"=============== \n";
	for(int i = 0; i <= (temp_p -> n); ++i){

		for(int j = 0; j <= (temp_p -> m); ++j){

			long long int index = Get_single_index(i, j, temp_p -> m + 1);

			std::cout<< i << " " << j <<" " <<temp_p -> solution[1][index] 
				<< " "<< temp_c -> solution[1][index]<< "\n" ;

		}
	}

}


#ifndef DG_UNIT_H
#define DG_UNIT_H

#include <vector>


/// @brief
/// Element unit. All the needed information is stored as a unit.
/// @param child_position ith-child
class Unit{

public:

	int n, m; 	// polynomial orders (x, y)	

	int index[3]{0, 0, 0}; 	// element index[i, j, k]
 
//	int local_index{};	// element local index	!! redundent?

	char status; 	// Hilbert status

	int child_position{};	// ith-child
	
	struct Face;	// forward declare
	std::vector<std::vector<Face>> facen;	// a structure for the neighbour on each face, initialized in dg_unit.cpp

	double xcoords[2]{0.0, 0.0};	// x coordinates
	double ycoords[2]{0.0, 0.0};	// y coordinates

	double* solution = nullptr;	// dymanic solution 

	Unit* next = nullptr;	// pointer to the next Unit
	

	// for testing-----------------------------------------
	int var{};
	double n_interface{};
	double s_interface{};
	//-----------------------------------------------------


	// constructor (default)
	Unit();
	
	// key function
//	int GetKey();
//
//	int Cantor_pairing(int x, int y);
};

/// @param face_type interface type of current direction.
/// @param hlevel hlevel of neighbour element.
/// @param porderx polynomial order of neighbour element in x direction.
/// @param pordery polynomial order of neighbour element in y direction.
/// @param key neighbour's key.
/// @param rank if face_type == 'M' rank == facing process's MPI.
struct Unit::Face{
	
	char face_type;
	
	int hlevel;
	
	int porderx;

	int pordery;
	
	int key;
	
	int rank;
};


#endif

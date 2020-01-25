#ifndef DG_UNIT_H
#define DG_UNIT_H

#include <vector>


/// @brief
/// Element unit. All the needed information is stored as a unit.
/// @param child_position ith-child
class Unit{

public:
	int level;	// h-refinement level !!Redundent? since we have index!!

	int n, m; 	// polynomial orders (x, y)	

	int index[3]{0, 0, 0}; 	// element index[i, j, k]
 
	int local_index{};	// element local index	!! redundent?

	char status; 	// Hilbert status

	int child_position{};	// ith-child
	
	short faces[4]{0, 0, 0, 0};	// remote element number (positive). Negtive means its on the physical boundary. 

	struct Face;	// forward declare
	std::vector<std::vector<Face>> facen;	// a structure for the neighbour on each face

	double xcoords[2]{0.0, 0.0};	// x coordinates
	double ycoords[2]{0.0, 0.0};	// y coordinates

	double* solution = nullptr;	// dymanic solution 

	Unit* next = nullptr;	// pointer to the next Unit
	

	// for testing-----------------------------------------
	int var{};
	//-----------------------------------------------------


	// constructor (default)
	Unit();
	
	// key function
//	int GetKey();
//
//	int Cantor_pairing(int x, int y);
};

struct Unit::Face{
	
	char face_type;
	
	int hlevel{};
	
	int porder{};
	
	int key{};
	
	int rank{};
};


#endif

#ifndef DG_UNIT_H
#define DG_UNIT_H

#include <vector>
#include <unordered_map>


/// @brief
/// Element unit. All the needed information is stored as a unit.
/// @param child_position ith-child
class Unit{

public:

	int n, m; 	// polynomial orders (x, y)	

	int index[3]{0, 0, 0}; 	// element index[i, j, k]
 
	char status; 	// Hilbert status

	int child_position{};	// ith-child
	
	struct Face;	// forward declare
	std::vector<std::vector<Face>> facen;	// a structure for the neighbour on each face, initialized in dg_unit.cpp

	double xcoords[2]{0.0, 0.0};	// x coordinates
	double ycoords[2]{0.0, 0.0};	// y coordinates

	std::unordered_map<int, std::vector<double>> solution;	// solutions

	std::vector<double> solution_int_l;	// solution on the element left interface
	std::vector<double> solution_int_r;	// solution on the element right interface

	std::unordered_map<int, std::vector<double>> ghost; // ghost space to store the neighbours' solutions on the interface
							    // hash <neighbour's key, solution_int>

	std::vector<double> nflux_l;	// numerical flux on the left interface
	std::vector<double> nflux_r;	// numerical flux on the right interface

	std::unordered_map<int, std::vector<double>> solution_time_der;	// solution time derivative

	Unit* next = nullptr;	// pointer to the next Unit

	bool hrefine = false;	// refinemnt 
	bool coarsen = false; 

	// for testing-----------------------------------------
	int var{};
	double n_interface{};
	double s_interface{};
	//-----------------------------------------------------


	// constructor (default)
	Unit();

		
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

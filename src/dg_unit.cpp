#include "dg_unit.h"
#include "dg_param.h"
#include <vector>


/// @brief
/// constructor (default).
/// default: level = 0, n = min poly order, m = max poly order
/// index = {0, 0, 0}, mpi_f = false, x/ycoords = {0.0, 0.0}
/// solution = nullptr. 
Unit::Unit() : n(grid::nmin), m(grid::nmin)
{
	facen = std::vector<std::vector<Face>>(4);

	ref_x = std::vector<double>(2);
	ref_y = std::vector<double>(2);

	// initialize element reference boudaries. Reference space [-1, 1]
	ref_x[0] = -1.0; ref_y[0] = -1.0;
	ref_x[1] =  1.0; ref_y[1] =  1.0;

}


// default construtor1
Unit::Face::Face() : hlevel{}, porderx{}, pordery{}, key{}, rank{}
{

	ref_x = std::vector<double> {-1.0, 1.0};
	ref_y = std::vector<double> {-1.0, 1.0};

}

// constructor2
Unit::Face::Face(char c, int h, int nx, int ny, int k, int r, double* ref1, double* ref2)
	: face_type(c), hlevel(h), porderx(nx), pordery(ny), key(k), rank(r)
{

	ref_x = std::vector<double> {ref1[0], ref1[1]};
	ref_y = std::vector<double> {ref2[0], ref2[1]};
}

// constructor3
Unit::Face::Face(char c, int h, int nx, int ny, int k, int r, std::vector<double>& ref1, std::vector<double>& ref2)
	: face_type(c), hlevel(h), porderx(nx), pordery(ny), key(k), rank(r)
{

	ref_x = ref1;
	ref_y = ref2;
}

// copy constructor ------------------------------------------------
Unit::Face::Face(const Face& face){	// copy another instance

	face_type = face.face_type;

	hlevel = face.hlevel;

	porderx = face.porderx;
	pordery = face.pordery;

	key = face.key;

	rank = face.rank;

	ref_x = face.ref_x;
	ref_y = face.ref_y;
	
}


Unit::Face::Face(const std::vector<Face>::iterator p){	// copy by pointer

	face_type = p -> face_type;

	hlevel = p -> hlevel;

	porderx = p -> porderx;
	pordery = p -> pordery;

	key = p -> key;

	rank = p -> rank;

	ref_x = p -> ref_x;
	ref_y = p -> ref_y;

}

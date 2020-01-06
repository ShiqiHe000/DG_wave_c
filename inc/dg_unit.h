#ifndef DG_UNIT_H
#define DG_UNIT_H

class Unit{

public:
	int level;	// h-refinement level

	int n, m; 	// polynomial orders (x, y)	

	int index[3]{0, 0, 0}; 	// element index[i, j, k]

	char status; 	// Hilbert status
	
	bool mpi_f[4]{false, false, false, false};	// mpi-boundary flag (true: on the mpi boundary)

	double xcoords[2]{0.0, 0.0};	// x coordinates
	double ycoords[2]{0.0, 0.0};	// y coordinates

	double* solution = nullptr;	// dymanic solution 
	
	
	// constructor (default)
	Unit();

};



#endif

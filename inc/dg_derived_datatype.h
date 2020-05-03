#ifndef DG_DERIVED_DATA_TYPE_H
#define DG_DERIVED_DATA_TYPE_H

#include <mpi.h>
#include <vector>

/// @brief
/// Packing the facen info and send it to neighbours to update MPI boundaries.
/// Use after hp-refinement. 
struct facen_pack{

	int local_key;	// neighbour element's key

	int hlevel;

	int porderx;
	int pordery;

	int mpi_length;

	double ref_x[2];
	double ref_y[2];

	// member function
	void Copy_ref(std::vector<double>& rx, std::vector<double>& ry);
};


/// @brief
/// Pack all the element information together and send together.
/// The variables inside the sturct refer to the same as class Unit.
struct info_pack{
	
	int n;	// assume n == m

	int index[3];

	char status;

	int child_position;

	double xcoords[2];
	double ycoords[2];

	double ref_x[2];
	double ref_y[2];
};

/// @brief
/// Pack all the element faces information together and send. 
/// Variables without annotation refer to the same one as class Unit. 
/// @param owners_key The key of the owner element.
/// @param facei face direction (0, 1, 2, 3). 
struct face_pack{

	int owners_key;	// the key of the owner element

	int facei; 	// face direction (0, 1, 2, 3)
	
	char face_type;
		
	int hlevel;
	
	int porderx;

	int pordery; 

	int key;

	int rank;

	double ref_x[2];
	double ref_y[2];

	// member function
	void Copy_ref(std::vector<double>& rx, std::vector<double>& ry);
};

namespace Hash{

	extern MPI_Datatype Facen_type;

	extern MPI_Datatype Elem_type;
	
	extern MPI_Datatype Face_type;
};


void Construct_data_type();
 
void Free_type();

#endif

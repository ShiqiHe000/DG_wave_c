#ifndef DG_DERIVED_DATA_TYPE_H
#define DG_DERIVED_DATA_TYPE_H

#include <mpi.h>
#include <stdio.h>

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

	int var;

};

/// @brief
/// Pack all the element faces information together and send. 
/// Variables without annotation refer to the same one as class Unit. 
/// @param owners_key The key of the owner element.
/// @param facei face direction (0, 1, 2, 3). 
//#pragma pack(push, 1)
struct face_pack{

	int owners_key;	// the key of the owner element

	int facei; 	// face direction (0, 1, 2, 3)
	
	char face_type;
		
	int hlevel;
	
	int porderx;

	int pordery; 

	int key;

	int rank;

};
//#pragma pack(pop)

namespace Hash{

	extern MPI_Datatype Elem_type;
	
	extern MPI_Datatype Face_type;
};


void Construct_data_type();
 
void Free_type();

#endif

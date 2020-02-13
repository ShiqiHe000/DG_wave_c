#ifndef DG_DERIVED_DATA_TYPE_H
#define DG_DERIVED_DATA_TYPE_H

#include <mpi.h>

/// @brief
/// Pack all the elment information together and send together.
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

namespace Hash{

	extern MPI_Datatype Elem_type;
};


void Construct_data_type();
 
void Free_type();

#endif

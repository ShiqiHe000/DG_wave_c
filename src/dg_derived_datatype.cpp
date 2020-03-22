#include <mpi.h>
#include "dg_derived_datatype.h"

namespace Hash{

	MPI_Datatype Elem_type;

	MPI_Datatype Face_type;
};

// forward declaration-----------------------------------------
void MPI_Elem_type();
void MPI_Face_type();
//--------------------------------------------------------------

/// @brief 
/// Construct data type to send the target element together.
void Construct_data_type(){

	MPI_Elem_type();

	MPI_Face_type();

}


/// @brief
/// Construct Face_type for sending info of element neighbours.  
void MPI_Face_type(){

	int num = 8;

	// Number of elements in each block (array of integers)
	int elem_blocklength[num]{1, 1, 1, 1, 1, 1, 1, 1};
	
	// Byte displacement of each block (array of integers).
	MPI_Aint array_of_offsets[num];
	MPI_Aint intex, charex;
	MPI_Aint lb;
	MPI_Type_get_extent(MPI_INT, &lb, &intex);
	MPI_Type_get_extent(MPI_CHAR, &lb, &charex);

	array_of_offsets[0] = (MPI_Aint) 0;	// owners_key
	array_of_offsets[1] = array_of_offsets[0] + intex ;	// facei
	array_of_offsets[2] = array_of_offsets[1] + intex ;	// face_type
	array_of_offsets[3] = array_of_offsets[2] + charex;	// hlevel
	array_of_offsets[4] = array_of_offsets[3] + intex;	// porderx
	array_of_offsets[5] = array_of_offsets[4] + intex;	// pordery
	array_of_offsets[6] = array_of_offsets[5] + intex;	// key
	array_of_offsets[7] = array_of_offsets[6] + intex;	// rank

	MPI_Datatype array_of_types[num]{MPI_INT, MPI_INT, MPI_CHAR, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};

	// create and MPI datatype
	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Face_type);	
	MPI_Type_commit(&Hash::Face_type);

	

}

/// @brief
/// Construct Elem_type for sending info of element variables.  
void MPI_Elem_type(){


	int num = 7;	// number of primitive MPI datatype

	// Number of elements in each block (array of integers)
	int elem_blocklength[num]{1, 3, 1, 1, 2, 2, 1};
	
	// Byte displacement of each block (array of integers).
	MPI_Aint array_of_offsets[num];
	MPI_Aint intex, charex, doubleex;
	MPI_Aint lb;
	MPI_Type_get_extent(MPI_INT, &lb, &intex);
	MPI_Type_get_extent(MPI_CHAR, &lb, &charex);
	MPI_Type_get_extent(MPI_DOUBLE, &lb, &doubleex);

	array_of_offsets[0] = (MPI_Aint) 0;
	array_of_offsets[1] = array_of_offsets[0] + intex;
	array_of_offsets[2] = array_of_offsets[1] + intex * 3;
	array_of_offsets[3] = array_of_offsets[2] + charex;
	array_of_offsets[4] = array_of_offsets[3] + intex;
	array_of_offsets[5] = array_of_offsets[4] + doubleex * 2;
	array_of_offsets[6] = array_of_offsets[5] + doubleex * 2;

	MPI_Datatype array_of_types[num]{MPI_INT, MPI_INT, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};

	// create and MPI datatype
	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Elem_type);	
	MPI_Type_commit(&Hash::Elem_type);


}

/// @brief
/// Free up the derived data type. 
void Free_type(){
	MPI_Type_free(&Hash::Elem_type);	

	MPI_Type_free(&Hash::Face_type);	

}

#include <mpi.h>
#include "dg_derived_datatype.h"
#include <vector>
#include <iostream>

namespace Hash{

	MPI_Datatype Elem_type;

	MPI_Datatype Face_type;
};

// forward declaration-----------------------------------------
void MPI_Elem_type();
void MPI_Face_type();
void MPI_Elem_type_new();
//--------------------------------------------------------------

/// @brief 
/// Construct data type to send the target element together.
void Construct_data_type(){

//	MPI_Elem_type();
	MPI_Elem_type_new();

	MPI_Face_type();

}

/// @brief
/// Construct Face_type for sending info of element neighbours.  
void MPI_Face_type(){

	int num = 3;

	int elem_blocklength[num]{2, 1, 5};

	MPI_Datatype array_of_types[num]{MPI_INT, MPI_CHAR, MPI_INT};
	
	MPI_Aint array_of_offsets[num];
	MPI_Aint baseadd, add1, add2;
	
	std::vector<face_pack> myface(1);

	MPI_Get_address(&(myface[0].owners_key), &baseadd);
	MPI_Get_address(&(myface[0].face_type), &add1);
	MPI_Get_address(&(myface[0].hlevel), &add2);

	array_of_offsets[0] = 0;
	array_of_offsets[1] = add1 - baseadd;
	array_of_offsets[2] = add2 - baseadd;

	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Face_type);	

	// check that the extent is correct
	MPI_Aint lb, extent;
	MPI_Type_get_extent(Hash::Face_type, &lb, &extent);	
	if(extent != sizeof(myface[0])){
		MPI_Datatype old = Hash::Face_type;
		MPI_Type_create_resized(old, 0, sizeof(myface[0]), &Hash::Face_type);
		MPI_Type_free(&old);
	}
	MPI_Type_commit(&Hash::Face_type);
}

void MPI_Elem_type_new(){

	int num = 5;

	int elem_blocklength[num]{4, 1, 1, 4, 1};
	
	MPI_Aint array_of_offsets[num];
	MPI_Aint baseadd, add1, add2, add3, add4;
	
	std::vector<info_pack> myinfo(1);

	MPI_Get_address(&(myinfo[0].n), &baseadd);
	MPI_Get_address(&(myinfo[0].status), &add1);
	MPI_Get_address(&(myinfo[0].child_position), &add2);
	MPI_Get_address(&(myinfo[0].xcoords), &add3);
	MPI_Get_address(&(myinfo[0].var), &add4);

	array_of_offsets[0] = 0;
	array_of_offsets[1] = add1 - baseadd;
	array_of_offsets[2] = add2 - baseadd;
	array_of_offsets[3] = add3 - baseadd;
	array_of_offsets[4] = add4 - baseadd;

	MPI_Datatype array_of_types[num]{MPI_INT, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_INT};

	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Elem_type);	

	MPI_Aint lb, extent;
	MPI_Type_get_extent(Hash::Elem_type, &lb, &extent);	
	if(extent != sizeof(myinfo[0])){
		MPI_Datatype old = Hash::Elem_type;
		MPI_Type_create_resized(old, 0, sizeof(myinfo[0]), &Hash::Elem_type);
		MPI_Type_free(&old);
	}
	MPI_Type_commit(&Hash::Elem_type);
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

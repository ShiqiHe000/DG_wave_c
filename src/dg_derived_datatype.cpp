#include <mpi.h>
#include "dg_derived_datatype.h"
#include <vector>
#include <iostream>

// memeber function of struct facen_pack ===========================================
void facen_pack::Copy_ref(std::vector<double>& rx, std::vector<double>& ry){
		
	for(int i = 0; i < 2; ++i){

		ref_x[i] = rx[i];
		ref_y[i] = ry[i];
	}
}
//===================================================================================

// memeber function of struct face_pack ===========================================
void face_pack::Copy_ref(std::vector<double>& rx, std::vector<double>& ry){
		
	for(int i = 0; i < 2; ++i){

		ref_x[i] = rx[i];
		ref_y[i] = ry[i];
	}
}
//===================================================================================

namespace Hash{

	MPI_Datatype Facen_type;

	MPI_Datatype Elem_type;

	MPI_Datatype Face_type;

	MPI_Datatype Adj_pairs;

	MPI_Datatype Owner_type;

};

// forward declaration-----------------------------------------
void MPI_Elem_type();
void MPI_Face_type();
void MPI_Facen_type();
void MPI_Adj_pairs_type();
void MPI_Owner_type();
//--------------------------------------------------------------

/// @brief 
/// Construct data type to send the target element together.
void Construct_data_type(){

	MPI_Facen_type();

	MPI_Elem_type();

	MPI_Face_type();
	
	MPI_Adj_pairs_type();

	MPI_Owner_type();

}


/// @brief
/// Construct a customized data type for exchange element ownership on the mpi boundaries. 
/// <long long int, int>
void MPI_Owner_type(){

	int num = 2;

	int elem_blocklength[num]{1, 1};

	MPI_Datatype array_of_types[num]{MPI_LONG_LONG_INT, MPI_INT};
	
	MPI_Aint array_of_offsets[num];
	MPI_Aint baseadd, add1;
	
	std::vector<owner_struct> my_owner(1);

	MPI_Get_address(&(my_owner[0].local_key), &baseadd);
	MPI_Get_address(&(my_owner[0].owners_rank), &add1);

	array_of_offsets[0] = 0;
	array_of_offsets[1] = add1 - baseadd;

	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Owner_type);	

	// check that the extent is correct
	MPI_Aint lb, extent;
	MPI_Type_get_extent(Hash::Owner_type, &lb, &extent);	
	if(extent != sizeof(my_owner[0])){
		MPI_Datatype old = Hash::Owner_type;
		MPI_Type_create_resized(old, 0, sizeof(my_owner[0]), &Hash::Owner_type);
		MPI_Type_free(&old);
	}
	MPI_Type_commit(&Hash::Owner_type);


}

/// @brief
/// Construct neighbour pairs.  
void MPI_Adj_pairs_type(){

	int num = 1;

	int elem_blocklength[num]{2};

	MPI_Datatype array_of_types[num]{MPI_LONG_LONG_INT};
	
	MPI_Aint array_of_offsets[num];
	MPI_Aint baseadd;
	
	std::vector<neighbour_pair> my_pair(1);

	array_of_offsets[0] = 0;

	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Adj_pairs);	

	// check that the extent is correct
	MPI_Aint lb, extent;
	MPI_Type_get_extent(Hash::Adj_pairs, &lb, &extent);	
	if(extent != sizeof(my_pair[0])){
		MPI_Datatype old = Hash::Adj_pairs;
		MPI_Type_create_resized(old, 0, sizeof(my_pair[0]), &Hash::Adj_pairs);
		MPI_Type_free(&old);
	}
	MPI_Type_commit(&Hash::Adj_pairs);


}


/// @brief
/// Construct Face_type for sending info of element neighbours.  
void MPI_Facen_type(){

	int num = 3;

	int elem_blocklength[num]{1, 4, 4};

	MPI_Datatype array_of_types[num]{MPI_LONG_LONG_INT, MPI_INT, MPI_DOUBLE};
	
	MPI_Aint array_of_offsets[num];
	MPI_Aint baseadd, add1, add2;
	
	std::vector<facen_pack> myface(1);

	MPI_Get_address(&(myface[0].local_key), &baseadd);
	MPI_Get_address(&(myface[0].hlevel), &add1);
	MPI_Get_address(&(myface[0].ref_x[0]), &add2);

	array_of_offsets[0] = 0;
	array_of_offsets[1] = add1 - baseadd;
	array_of_offsets[2] = add2 - baseadd;

	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Facen_type);	

	// check that the extent is correct
	MPI_Aint lb, extent;
	MPI_Type_get_extent(Hash::Facen_type, &lb, &extent);	
	if(extent != sizeof(myface[0])){
//std::cout << "extent " << extent << " size " << sizeof(myface[0]) << "\n";
		MPI_Datatype old = Hash::Facen_type;
		MPI_Type_create_resized(old, 0, sizeof(myface[0]), &Hash::Facen_type);
		MPI_Type_free(&old);
	}
	MPI_Type_commit(&Hash::Facen_type);
}



/// @brief
/// Construct Face_type for sending info of element neighbours.  
void MPI_Face_type(){

	int num = 7;

	int elem_blocklength[num]{1, 1, 1, 3, 1, 1, 4};

	MPI_Datatype array_of_types[num]{MPI_LONG_LONG_INT, MPI_INT, MPI_CHAR, MPI_INT, MPI_LONG_LONG_INT, MPI_INT, MPI_DOUBLE};
	
	MPI_Aint array_of_offsets[num];
	MPI_Aint baseadd, add1, add2, add3, add4, add5, add6;
	
	std::vector<face_pack> myface(1);

	MPI_Get_address(&(myface[0].owners_key), &baseadd);
	MPI_Get_address(&(myface[0].facei), &add1);
	MPI_Get_address(&(myface[0].face_type), &add2);
	MPI_Get_address(&(myface[0].hlevel), &add3);
	MPI_Get_address(&(myface[0].key), &add4);
	MPI_Get_address(&(myface[0].rank), &add5);
	MPI_Get_address(&(myface[0].ref_x[0]), &add6);

	array_of_offsets[0] = 0;
	array_of_offsets[1] = add1 - baseadd;
	array_of_offsets[2] = add2 - baseadd;
	array_of_offsets[3] = add3 - baseadd;
	array_of_offsets[4] = add4 - baseadd;
	array_of_offsets[5] = add5 - baseadd;
	array_of_offsets[6] = add6 - baseadd;

	MPI_Type_create_struct(num, elem_blocklength, array_of_offsets, array_of_types, &Hash::Face_type);	

	// check that the extent is correct
	MPI_Aint lb, extent;
	MPI_Type_get_extent(Hash::Face_type, &lb, &extent);	
	if(extent != sizeof(myface[0])){
//std::cout << "=====================prob \n";
		MPI_Datatype old = Hash::Face_type;
		MPI_Type_create_resized(old, 0, sizeof(myface[0]), &Hash::Face_type);
		MPI_Type_free(&old);
	}
	MPI_Type_commit(&Hash::Face_type);
}

void MPI_Elem_type(){

	int num = 4;

	int elem_blocklength[num]{4, 1, 1, 8};
	
	MPI_Aint array_of_offsets[num];
	MPI_Aint baseadd, add1, add2, add3; 
	
	std::vector<info_pack> myinfo(1);

	MPI_Get_address(&(myinfo[0].n), &baseadd);
	MPI_Get_address(&(myinfo[0].status), &add1);
	MPI_Get_address(&(myinfo[0].child_position), &add2);
	MPI_Get_address(&(myinfo[0].xcoords[0]), &add3);

	array_of_offsets[0] = 0;
	array_of_offsets[1] = add1 - baseadd;
	array_of_offsets[2] = add2 - baseadd;
	array_of_offsets[3] = add3 - baseadd;

	MPI_Datatype array_of_types[num]{MPI_INT, MPI_CHAR, MPI_INT, MPI_DOUBLE};

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
/// Free up the derived data type. 
void Free_type(){
	MPI_Type_free(&Hash::Facen_type);	

	MPI_Type_free(&Hash::Elem_type);	

	MPI_Type_free(&Hash::Face_type);	

	MPI_Type_free(&Hash::Adj_pairs);	

	MPI_Type_free(&Hash::Owner_type);	

}

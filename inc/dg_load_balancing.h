#ifndef DG_LOAD_BALANCING_H
#define DG_LOAD_BALANCING_H

/// @brief
/// This structure assist to form the processor mapping table.
/// @param iproc Process's rank number.
/// @param gnum first element's global number.
/// @param key first element's key.
struct pmap{

	int irank;	

	int gnum;

//	int key;

};

void Build_mapping_table();


#endif

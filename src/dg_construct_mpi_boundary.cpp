#include "dg_construct_mpi_boundary.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_cantor_pairing.h"
#include "dg_nodal_2d_storage.h"
#include <unordered_map>

/// @brief
/// Construct MPI boundaries and physical boundaries. 
/// Each element has 4 faces
/// \verbatim
///   ---f2----
///   |       |
/// f3|       |f4
///   |       |
///   ---f1----
/// \endverbatim
/// Boundaries circumstances are stored in short face[4]. 
/// positive value means the number of element on the other side of the MPI boundary.
/// negative value means the face in on the physical boundary and the magnitude of 
/// the value represents the face number as shown above. 
void MPI_boundary_construct(){
	
	// start by the first element
	Unit* temp = local::head;
	
	// traverse the linked-list
	for(int k = 0; k < local::local_elem_num; ++k){

		// on the south physical boundary?
		if(temp -> index[0] == 0){
			temp -> faces[0] = -1;	// Yes

			temp -> facen[0].push_back(Unit::Face());

			temp -> facen[0][0].face_type = 'B'; // Yes
		}
		else{	// No, search south neighbour
			int ni = temp -> index[0] - 1;
			int nj = temp -> index[1];
			int nkey = Get_key_fun(ni, nj, 0); 	// before adapt

			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[0] = 1;
			}
		}

		// on the north physical boundary?
		if(temp -> index[0] == (SortMesh::num_of_element_x - 1)){
			temp -> faces[1] = -2;	// yes

		}
		else{	// no
			int ni = temp -> index[0] + 1;
			int nj = temp -> index[1];
			int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[0] = 1;
			}
			
		}

		// on the west physical boundary?
		if(temp -> index[1] == 0){
			temp -> faces[2] = -3;	// yes

		}
		else{	// no
			int ni = temp -> index[0];
			int nj = temp -> index[1] - 1;
			int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[2] = 1;
			}
			
		}

		// on the east physical boundary?
		if(temp -> index[1] == SortMesh::num_of_element_y - 1){
			temp -> faces[3] = -4;	// yes

		}
		else{	// no
			int ni = temp -> index[0];
			int nj = temp -> index[1] + 1;
			int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[3] = 1;
			}
			
		}
		
		// pointer move to next element
		temp = temp -> next;
	}

}

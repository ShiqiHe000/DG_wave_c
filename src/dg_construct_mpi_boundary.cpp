#include "dg_construct_mpi_boundary.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_cantor_pairing.h"
#include "dg_nodal_2d_storage.h"
#include <unordered_map>
#include "dg_param.h"
#include "dg_search_rank.h"

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

			temp -> facen[0].push_back(Unit::Face());

			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[0] = 1;
			
				temp -> facen[0][0].face_type = 'M';
				temp -> facen[0][0].porder = grid::nmin;	// for uniform mesh, we can record
				temp -> facen[0][0].key = nkey;		// hlevel initially 0

				int target_rank = Target_rank(ni, nj);

				temp -> facen[0][0].rank = target_rank;
				
			}
			else{	// if found, record info

				temp -> facen[0][0].face_type = 'L';
				temp -> facen[0][0].porder = grid::nmin;	
				temp -> facen[0][0].key = nkey;	
			

			}
		}

		// on the north physical boundary?
		temp -> facen[1].push_back(Unit::Face());
		if(temp -> index[0] == (SortMesh::num_of_element_x - 1)){
			temp -> faces[1] = -2;	// yes
			
			temp -> facen[1][0].face_type = 'B';

		}
		else{	// no
			int ni = temp -> index[0] + 1;
			int nj = temp -> index[1];
			int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[1] = 1;

				temp -> facen[1][0].face_type = 'M';
				temp -> facen[1][0].porder = grid::nmin;	
				temp -> facen[1][0].key = nkey;	
				int target_rank = Target_rank(ni, nj);

				temp -> facen[1][0].rank = target_rank;
			}
			else{

				temp -> facen[1][0].face_type = 'L';
				temp -> facen[1][0].porder = grid::nmin;	
				temp -> facen[1][0].key = nkey;	

			}
			
		}

		// on the west physical boundary?
		temp -> facen[2].push_back(Unit::Face());
		if(temp -> index[1] == 0){
			temp -> faces[2] = -3;	// yes

			temp -> facen[2][0].face_type = 'B';

		}
		else{	// no
			int ni = temp -> index[0];
			int nj = temp -> index[1] - 1;
			int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[2] = 1;
				
				temp -> facen[2][0].face_type = 'M';
				temp -> facen[2][0].porder = grid::nmin;	
				temp -> facen[2][0].key = nkey;	
				int target_rank = Target_rank(ni, nj);

				temp -> facen[2][0].rank = target_rank;
			}
			else{
				
				temp -> facen[2][0].face_type = 'L';
				temp -> facen[2][0].porder = grid::nmin;	
				temp -> facen[2][0].key = nkey;	

			}
			
		}

		// on the east physical boundary?
		temp -> facen[3].push_back(Unit::Face());
		if(temp -> index[1] == SortMesh::num_of_element_y - 1){
			temp -> faces[3] = -4;	// yes

			temp -> facen[3][0].face_type = 'B';
		}
		else{	// no
			int ni = temp -> index[0];
			int nj = temp -> index[1] + 1;
			int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				temp -> faces[3] = 1;
				
				temp -> facen[3][0].face_type = 'M';
				temp -> facen[3][0].porder = grid::nmin;	
				temp -> facen[3][0].key = nkey;	

				int target_rank = Target_rank(ni, nj);

				temp -> facen[3][0].rank = target_rank;
			}
			else{

				temp -> facen[3][0].face_type = 'L';
				temp -> facen[3][0].porder = grid::nmin;	
				temp -> facen[3][0].key = nkey;	
				
			
			}
			
		}
		
		// pointer move to next element
		temp = temp -> next;
	}

}

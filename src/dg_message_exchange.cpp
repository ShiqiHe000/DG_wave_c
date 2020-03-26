#include <mpi.h>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include <vector>
#include <unordered_map>

/// @brief
/// Exchange element interface info for those on the MPI boundaries. 
/// North send, south recv.
/// @param
void Exchange_solution(std::unordered_map<int, std::vector<int>>& sender, std::unordered_map<int, std::vector<int>>& recver){



}

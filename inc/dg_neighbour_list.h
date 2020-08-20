#ifndef DG_NEIGHBOUR_LIST
#define DG_NEIGHBOUR_LIST

#include <unordered_map>
#include <vector>

void Neighbours_array_x(int i, int j, int k, long long int local_key, int facen, 
			std::unordered_map<long long int, std::vector<long long int>>& neighbours);

void Neighbours_array_y(int i, int j, int k, long long int local_key, int facen, 
			std::unordered_map<long long int, std::vector<long long int>>& neighbours);
#endif

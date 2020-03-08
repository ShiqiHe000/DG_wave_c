#ifndef DG_NEIGHBOUR_LIST
#define DG_NEIGHBOUR_LIST

void Total_neighbours(int level_now, int& all, int& my_position);

void Neighbours_array_x(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours);

void Neighbours_array_y(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours);
#endif

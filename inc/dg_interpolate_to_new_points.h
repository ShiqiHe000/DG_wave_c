#ifndef DG_INTERPOLATE_TO_NEW_POINT_H
#define DG_INTERPOLATE_TO_NEW_POINT_H

#include <array>
#include <vector>

void Polynomial_interpolate_matrix(std::vector<double>& x, std::vector<double>& xi, std::vector<double>& T);

void Solutions_to_children(std::array<long long int, 4>& keys, long long int p_key);

void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start_old, int start_new, int interval);

void Solution_back_to_parent(std::array<long long int, 4>& keys, long long int p_key);

#endif

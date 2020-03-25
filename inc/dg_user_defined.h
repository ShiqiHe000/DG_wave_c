#ifndef DG_USER_DEFINED_H
#define DG_USER_DEFINED_H

#include <vector>
#include <unordered_map>

void Exact_solution_Gaussian(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t);
#endif

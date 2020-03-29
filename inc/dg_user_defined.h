#ifndef DG_USER_DEFINED_H
#define DG_USER_DEFINED_H

#include <vector>
#include <unordered_map>

// user defined variable-------------------------------
namespace user{
	extern const double kx; 
	extern const double ky; 
	extern const double D;
	extern const double xx0;
	extern const double yy0; 
};
//------------------------------------------------------

void Exact_solution_Gaussian(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t);
#endif

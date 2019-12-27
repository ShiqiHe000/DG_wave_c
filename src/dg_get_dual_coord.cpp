#include <mpi.h>
#include <cmath>

void Get_dual_coord_2d(double& x, double& y, double& delta_x, double& delta_y, int* coord){

	coord[0] = std::round(x / delta_x); 
	coord[1] = std::round(x / delta_y); 


}

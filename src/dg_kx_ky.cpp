#include <cmath>

void Get_kx_ky(double& kx, double& ky, double t){

	const double pi = 4.0 * std::atan(1.0);

	kx = sin(10000 * pi * t);
	ky = cos(10000 * pi * t);

}

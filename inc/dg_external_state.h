#ifndef DG_EXTERNAL_STATE_H
#define DG_EXTERNAL_STATE_H

#include <vector>

void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_sin_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_reflect_x(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);

void External_mirror_y_left(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_reflect_y(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);
#endif

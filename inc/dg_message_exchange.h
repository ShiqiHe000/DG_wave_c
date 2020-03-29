#ifndef DG_MESSAGE_EXCHANEG_H
#define DG_MESSAGE_EXCHANGE_H

#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"

void Exchange_solution(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s,
			std::unordered_map<int, std::vector<mpi_table>>& recver, int face_r, char dir);


void Exchange_flux(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s, int face_r);
#endif

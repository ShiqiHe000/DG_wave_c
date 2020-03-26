#ifndef DG_MPI_TABLE_CONSTRUCT_H
#define DG_MPI_TABLE_CONSTRUCT_H

#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"

void Construct_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& north, int face_north, 
				std::unordered_map<int, std::vector<int>>& neighbours_north, 
				std::unordered_map<int, std::vector<mpi_table>>& south, int face_south, 
				std::unordered_map<int, std::vector<int>>& neighbours_south);

void Update_mpi_boundaries(std::unordered_map<int, std::vector<mpi_table>>& north, int facen,
				std::unordered_map<int, std::vector<int>>& neighbours_north,
				std::unordered_map<int, std::vector<mpi_table>>& south, int faces, 
				std::unordered_map<int, std::vector<int>>& neighbours_south);

void Clear_tables();
#endif

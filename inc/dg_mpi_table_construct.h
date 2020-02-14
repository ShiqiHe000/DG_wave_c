#ifndef DG_MPI_TABLE_CONSTRUCT_H
#define DG_MPI_TABLE_CONSTRUCT_H

#include "dg_boundary_table.h"
#include <vector>

void Construct_mpi_table(std::vector<table_elem>& north, int facen, std::vector<table_elem>& south, int faces);

void Update_mpi_boundaries(std::vector<table_elem>& north, int facen, std::vector<accum_elem>& north_accum, 
				std::vector<table_elem>& south, int faces, std::vector<accum_elem>& south_accum);

void Clear_tables();
#endif

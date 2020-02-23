#ifndef DG_MPI_TABLE_CONSTRUCT_H
#define DG_MPI_TABLE_CONSTRUCT_H

#include "dg_boundary_table.h"
#include <vector>

void Construct_mpi_table_x(std::vector<table_elem>& north, std::vector<table_elem>& south);
void Construct_mpi_table_y(std::vector<table_elem>& west, std::vector<table_elem>& east);

void Update_mpi_boundaries(std::vector<table_elem>& north, int facen, std::vector<accum_elem>& north_accum, 
				std::vector<table_elem>& south, int faces, std::vector<accum_elem>& south_accum);

void Clear_tables();
#endif

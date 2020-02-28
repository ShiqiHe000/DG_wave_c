#ifndef DG_MPI_TABLE_CONSTRUCT_H
#define DG_MPI_TABLE_CONSTRUCT_H

#include "dg_boundary_table.h"
#include <vector>

void Construct_mpi_table_x();
void Construct_mpi_table_y();

void Update_mpi_boundaries(std::unordered_map& north, int facen, std::vector<table_elem>& south, int faces);

void Clear_tables();
#endif

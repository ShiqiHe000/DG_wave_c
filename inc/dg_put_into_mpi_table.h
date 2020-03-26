#ifndef DG_PUT_INTO_MPI_TABLE_H
#define DG_PUT_INTO_MPI_TABLE_H

#include "dg_boundary_table.h"

void Put_in_mpi_table(Unit* temp, std::vector<Unit::Face>::iterator& facen_it, 
			std::unordered_map<int, std::vector<mpi_table>>& table);

#endif


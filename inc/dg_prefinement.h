#ifndef DG_PREFINEMENT_H
#define DG_PREFINEMENT_H

#include "dg_unit.h"

void p_refinement_apply(Unit* temp);

void p_refinement(int kt);

void p_coarsening_interpolate(Unit* temp);

void p_coarsening_L2(Unit* temp);

void Update_neighbours_facen(Unit* temp);

#endif

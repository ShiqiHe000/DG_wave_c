#include "dg_error_estimator.h"
#include "dg_unit.h"
#include <vector>
#include <cassert>

/// @brief
/// 
/// @param num number of point to fit regression.
/// @param temp pointer to the estimated element.
double Error_indicator(int num, Unit* temp){

	assert((temp -> n > 6) && "Polynomial order is too low to estimate error.");

}

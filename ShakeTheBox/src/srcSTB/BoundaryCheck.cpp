/*
 * BoundaryCheck.cpp
 *
 *  Created on: May 16, 2018
 *      Author: Shiyong Tan
 */

#include "BoundaryCheck.h"

void BoundaryCheck::SetLimit(double x_ulim, double x_llim, double y_ulim, double y_llim,
		double z_ulim, double z_llim) {
	x_upper_lim = x_ulim; x_lower_lim = x_llim;
	y_upper_lim = y_ulim; y_lower_lim = y_llim;
	z_upper_lim = z_ulim; z_lower_lim = z_llim;
}

bool BoundaryCheck::Check(Position pos) {
	if (pos.X() <= x_upper_lim && pos.X() >= x_lower_lim
			&& pos.Y() <= y_upper_lim && pos.Y() >= y_lower_lim
			&& pos.Z() <= z_upper_lim && pos.Z() >= z_lower_lim) return true;
	return false;
}

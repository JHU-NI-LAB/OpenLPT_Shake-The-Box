/*
 * BoundaryCheck.h
 *
 *  Created on: May 16, 2018
 *      Author: Shiyong Tan
 */

#ifndef BOUNDARYCHECK_H_
#define BOUNDARYCHECK_H_

#include "Position.h"

class BoundaryCheck {
public:
	BoundaryCheck() {};
	BoundaryCheck(double x_ulim, double x_llim, double y_ulim, double y_llim,
			double z_ulim, double z_llim) : x_upper_lim(x_ulim), x_lower_lim(x_llim),
					y_upper_lim(y_ulim), y_lower_lim(y_llim), z_upper_lim(z_ulim), z_lower_lim(z_llim) {};

	~BoundaryCheck() {};

	void SetLimit(double x_ulim, double x_llim, double y_ulim, double y_llim,
			double z_ulim, double z_llim);

	bool Check(Position pos);

private:
	double x_upper_lim, x_lower_lim;
	double y_upper_lim, y_lower_lim;
	double z_upper_lim, z_lower_lim;
};


#endif /* BOUNDARYCHECK_H_ */

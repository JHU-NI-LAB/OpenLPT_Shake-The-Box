/*
 *  Camera.cpp
 *  
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 */

#include <Camera.h>
#include <string.h>
#include <limits.h>

using namespace std;

Camera::Camera(istream& is)
{
    is >> Noffh;
    is >> Noffw;
	is >> Npixw;
	is >> Npixh;
	is >> wpix;
	is >> hpix;
	is >> f_eff;
	is >> kr;
	is >> kx;
	double buffer[9];
	for (int i = 0; i < 9; ++i) {
		is >> buffer[i];
	}
	R = Matrix(buffer);
	for (int i = 0; i < 3; ++i) {
		is >> buffer[i];
	}
	T = Position(buffer[0], buffer[1], buffer[2]);
	for (int i = 0; i < 9; ++i) {
		is >> buffer[i];
	}
	Rinv = Matrix(buffer);
	for (int i = 0; i < 3; ++i) {
		is >> buffer[i];
	}
	Tinv = Position(buffer[0], buffer[1], buffer[2]);
}

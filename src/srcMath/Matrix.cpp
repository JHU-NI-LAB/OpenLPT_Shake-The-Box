/*
 *  Matrix.cpp
 *  
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 */

#include <Matrix.h>
#include <string.h>
#include <limits.h>

using namespace std;

Matrix Matrix::Invert()
{
	// this is messy!
	
	double determinant = (entries[0] * (entries[4]*entries[8] - entries[5]*entries[7])
												- entries[1] * (entries[3]*entries[8] - entries[5]*entries[6])
												+ entries[2] * (entries[3]*entries[7] - entries[4]*entries[6]));
	
	Matrix inverse;
	inverse.entries[0] = (entries[4]*entries[8] - entries[5]*entries[7]) / determinant;
	inverse.entries[1] = (entries[2]*entries[7] - entries[1]*entries[8]) / determinant;
	inverse.entries[2] = (entries[1]*entries[5] - entries[2]*entries[4]) / determinant;
	inverse.entries[3] = (entries[5]*entries[6] - entries[3]*entries[8]) / determinant;
	inverse.entries[4] = (entries[0]*entries[8] - entries[2]*entries[6]) / determinant;
	inverse.entries[5] = (entries[2]*entries[3] - entries[0]*entries[5]) / determinant;
	inverse.entries[6] = (entries[3]*entries[7] - entries[4]*entries[6]) / determinant;
	inverse.entries[7] = (entries[1]*entries[6] - entries[0]*entries[7]) / determinant;
	inverse.entries[8] = (entries[0]*entries[4] - entries[1]*entries[3]) / determinant;
	
	return inverse;
}

const Position operator*(const Position& p, const Matrix& m) 
{
	return Position(p.X() * m.entries[0] + p.Y() * m.entries[3] + p.Z() * m.entries[6],
									p.X() * m.entries[1] + p.Y() * m.entries[4] + p.Z() * m.entries[7],
									p.X() * m.entries[2] + p.Y() * m.entries[5] + p.Z() * m.entries[8]);
}

const Position operator*(const Matrix& m, const Position& p)
{
	return Position(m.entries[0] * p.X() + m.entries[1] * p.Y() + m.entries[2] * p.Z(),
									m.entries[3] * p.X() + m.entries[4] * p.Y() + m.entries[5] * p.Z(),
									m.entries[6] * p.X() + m.entries[7] * p.Y() + m.entries[8] * p.Z());
}

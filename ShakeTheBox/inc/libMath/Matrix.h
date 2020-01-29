/*
 *  Matrix.h
 *  
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>

#include <Position.h>
#include <string.h>
#include <limits.h>

class Matrix {
public:
	
	Matrix();
	Matrix(double buffer[]);
	Matrix(const Matrix& m);
	~Matrix() {};
	
	double Get(int i, int j) const;
	void Set(int i, int j, double val);
	
	Matrix Invert();
	
	Matrix& operator+=(const Matrix& right);
	
	friend const Position operator*(const Position& p, const Matrix& m);
	friend const Position operator*(const Matrix& m, const Position& p);
	
private:
	double entries[9];
};

inline Matrix::Matrix()
{
	for (int i = 0; i < 9; ++i) {
		entries[i] = 0;
	}
}

inline Matrix::Matrix(double buffer[])
{
	memcpy(entries, buffer, 9 * sizeof(double));
}

inline Matrix::Matrix(const Matrix& m)
{
	memcpy(entries, m.entries, 9 * sizeof(double));
}

inline double Matrix::Get(int i, int j) const
{
	return entries[i * 3 + j];
}

inline void Matrix::Set(int i, int j, double val)
{
	entries[i * 3 + j] = val;
}

inline Matrix& Matrix::operator+=(const Matrix& right)
{
	for (int i = 0; i < 9; ++i) {
		entries[i] += right.entries[i];
	}
	
	return *this;
}

#endif // MATRIX_H

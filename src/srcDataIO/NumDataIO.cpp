/*
 * NumDataIO.cpp
 *
 *  Created on: Feb 6, 2018
 *      Author: shiyongtan
 */
#include <fstream>
#include <sstream>
#include <string>
#include "NumDataIO.h"

/*
 * The code of template class in this file is insufficient for compiler to produce the methods
 * in this file which are needed in other file. So it would cause undefined reference problem.
 * there is two solutions:
 * 1. to put all the code in this file in the header file.
 * 2. to instantiate all possible templates as follows in the end of source file:
 * template class NumDataIO<int>;
 * template class NumDataIO<double>;
 */

// save data to txt file
template <class T>
int NumDataIO <T>::WriteData(T* data) {
	std::ofstream outfile;
	outfile.open(this->m_file_path);
	for (int i = 0; i < m_total_number; i++ ) {
		outfile<<*(data + i)<<","; //using comma as the delimiter.
	}
	outfile.close();
	return 0;
}

// read data from txt file
template <class T>
int NumDataIO <T>::ReadData(T* data) {
	std::ifstream infile;
	infile.open(this->m_file_path);
	std::string line;
	std::getline(infile, line); // reading a line into line.
	std::stringstream ss(line); // make the line as a stringstream to separate the line
	std::string cell; //to save each element
	for (int i = 0; i < m_total_number; i++ ) {
		getline(ss, cell, ',');
		*(data + i) = stod(cell);
	}
	infile.close();
	return 0;
}

template class NumDataIO<int>;
template class NumDataIO<double>;



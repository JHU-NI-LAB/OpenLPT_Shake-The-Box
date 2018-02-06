/*
 * DoubleDataIO.c
 *
 *  Created on: Feb 2, 2018
 *      Author: shiyongtan
 */

#include <fstream>
#include <sstream>
#include <string>
#include "DoubleDataIO.h"

// save data to txt file
int DoubleDataIO::WriteData(double* data) {
	std::ofstream outfile;
	outfile.open(m_file_path);
	for (int i = 0; i < m_total_number; i++ ) {
		outfile<<*(data + i)<<","; //using comma as the delimiter.
	}
	outfile.close();
	return 0;
}

// read data to txt file
int DoubleDataIO::ReadData(double* data) {
	std::ifstream infile;
	infile.open(m_file_path);
	std::string line;
	std::getline(infile, line); // reading a line into line.
	std::stringstream ss(line); // make the line as a stringstream to separate the line
	std::string cell; //to save each element
	for (int i = 0; i < m_total_number; i++ ) {
		getline(ss, cell, ',');
		*(data + i) = stod(cell);
	}
	return 0;
}






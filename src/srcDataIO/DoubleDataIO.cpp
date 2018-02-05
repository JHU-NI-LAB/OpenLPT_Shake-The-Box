/*
 * DoubleDataIO.c
 *
 *  Created on: Feb 2, 2018
 *      Author: shiyongtan
 */

#include <fstream>
#include "DoubleDataIO.h"

// save data to txt file
int DoubleDataIO::WriteData(double* data) {
	std::ofstream outfile;
	outfile.open(m_file_path);
	for (int i = 0; i < m_total_number; i++ ) {
		outfile<<*(data + i)<<"\n";
	}
	outfile.close();
	return 0;
}

// read data to txt file
int DoubleDataIO::ReadData(double* data) {
	std::ifstream infile;
	infile.open(m_file_path);
	for (int i = 0; i < m_total_number; i++ ) {
		infile>>*(data + i);
	}
	return 0;
}





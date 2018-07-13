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
#include "Common.h"

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
	if (m_save_mode == 1) {
		outfile.open(this->m_file_path, std::ios::app);
	} else {
		outfile.open(this->m_file_path);
	}
	outfile.precision(m_num_precision);
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
	if (infile.is_open()) {
		std::string line;
		std::getline(infile, line); // reading a line into line.
		std::stringstream ss(line); // make the line as a stringstream to separate the line
		std::string cell; //to save each element
		if (m_total_number <= 0) m_total_number = GetTotalNumber();
		for (unsigned int i = 0; i < m_total_number + m_skip_data_num; i++ ) {
			getline(ss, cell, ',');
			if (i >= m_skip_data_num) {  // skip the set number of data
				if (std::is_integral<T>::value) { // pass the value of cell according to different data type
					*(data + i - m_skip_data_num) = stoi(cell);
				} else {
					try {
						*(data + i - m_skip_data_num) = stod(cell);
					}
					catch (const std::out_of_range& oor) {
						*(data + i - m_skip_data_num) = 0;
					}

				}
			}
		}
		infile.close();
	} else {
		error = NO_FILE;
	}
	return 0;
}

template<class T>
int NumDataIO<T>::GetTotalNumber() {
	std::ifstream infile;
	infile.open(this->m_file_path);
	if (!infile.is_open()) {
		error = NO_FILE;
		return 0;
	}
	std::string line;
	std::getline(infile, line); // reading a line into line.
	std::stringstream ss(line); // make the line as a stringstream to separate the line
	std::string cell; //to save each element
	int total_number = 0;
	while (getline(ss, cell, ',')) total_number++;
	return total_number;
}

template class NumDataIO<int>;
template class NumDataIO<double>;



/*
 * DataIO.h
 *
 *  Created on: Feb 2, 2018
 *      Author: shiyongtan
 *  Description: Given txt file path, create or load data from/to the file.
 *
 */

#ifndef DATAIO_H_
#define DATAIO_H_

#include <string>

template <class T> class DataIO {
public:
	DataIO () {};
	~DataIO () {};

	/*
	 * Function: a virtual function to save data.
	 * 			its realization depends on different data type.
	 * Input: Data to be saved to the file.
	 * Output: An int to indicate the error state.
	 */
	virtual int WriteData(T* data) = 0;

	/*
	 * Function: a virtual function to read data.
	 * 			its realization depends on different data type.
	 * Input: Data to save the data from file.
	 * Output: An int to indicate the error state.
	 */
	virtual int ReadData(T* data) = 0;

	/*
	 * Function: set the file path
	 * Input: A string of file_path
	 * output: None
	 */
	void SetFilePath(std::string file_path)
	{ m_file_path = file_path; };

protected:
	std::string m_file_path;	// file path
};




#endif /* DATAIO_H_ */

/*
 * DoubleDataIO.h
 *
 *  Created on: Feb 2, 2018
 *      Author: shiyongtan
 *      Description: This derivative function is used to save the double type data only.
 *      ATTENTION: data should be defined clearly before using this class method.
 */

#ifndef NUMDATAIO_H_
#define NUMDATAIO_H_

#include "DataIO.h"

template <class T>
class NumDataIO : public DataIO<T> {
public:
	NumDataIO() {};
	~NumDataIO() {};

	/*
	 * Function: a realization to the virtual function
	 * 			 to save double type data to the txt file
	 * Input: Data to be saved to the file.
	 * Output: An int to indicate the error state.
	 */
	int WriteData(T* data);

	/*
	 * Function: a realization to the virtual function
	 * 			to read double type data from the txt file
	 * Input: Data to save the data from file.
	 * Output: An int to indicate the error state.
	 */
	int ReadData(T* data);

	/*
	 * Function: set the total number of elements.
	 * Input: An int of the number of elements.
	 * Output: None.
	 */
	void SetTotalNumber(int total_number) {
		m_total_number = total_number;
	}

	/*
	 * Function: get the total number of elements in txt file.
	 * Input: None
	 * Output: the total number.
	 */
	int GetTotalNumber();

	/*
	 * Function: to set saving mode, like to overwrite or to append
	 * Input: integer 0 means to overwrite, 1 to append
	 * Output: None
	 */
	void SaveMode(int mode) {
		m_save_mode = mode;
	}

	/*
	 * Function: to set the number of data to be skipped
	 * Input: int to indicate the number of data to be skipped
	 * Output: None.
	 */
	void SetSkipDataNum(int num) {
		m_skip_data_num = num;
	}

	/*
	 * Function: to set the precision of data
	 * Input: int for the precision
	 * Output: None.
	 */
	void SetNumPrecsion(int num) {
		m_num_precision = num;
	}

private:
	int m_total_number = 0;	// total number of the matrix elements
	int m_save_mode = 0;    // save mode:0 means to overwrite, 1 to append
	int m_skip_data_num = 0; // the number of data to be skipped
	int m_num_precision = 8; // the precision to save the data.
};




#endif /* NUMDATAIO_H_ */



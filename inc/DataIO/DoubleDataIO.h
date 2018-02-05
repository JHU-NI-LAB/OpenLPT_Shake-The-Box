/*
 * DoubleDataIO.h
 *
 *  Created on: Feb 2, 2018
 *      Author: shiyongtan
 *      Description: This derivative function is used to save the double type data only.
 *      ATTENTION: data should be defined clearly before using this class method.
 */

#ifndef DOUBLEDATAIO_H_
#define DOUBLEDATAIO_H_

#include "DataIO.h"

class DoubleDataIO : public DataIO<double> {
public:
	DoubleDataIO(){} ;
	~DoubleDataIO() {};

	/*
	 * Function: a realization to the virtual function
	 * 			 to save double type data to the txt file
	 * Input: Data to be saved to the file.
	 * Output: An int to indicate the error state.
	 */
	int WriteData(double* data);

	/*
	 * Function: a realization to the virtual function
	 * 			to read double type data from the txt file
	 * Input: Data to save the data from file.
	 * Output: An int to indicate the error state.
	 */
	int ReadData(double* data);

	/*
	 * Function: set the total number of elements.
	 * Input: An int of the number of elements.
	 * Output: None.
	 */
	void SetTotalNumber(int total_number) {
		m_total_number = total_number;
	}


private:
	int m_total_number;	// total number of the matrix elements
};



#endif /* DOUBLEDATAIO_H_ */

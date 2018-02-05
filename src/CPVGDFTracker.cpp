/*
 *
 */

#include "DoubleDataIO.h"
#include <string>

int main() {

	std::string file_path = "/home/shiyongtan/Documents/Code/TestFile/test.txt";
	DoubleDataIO double_data_io;
	double_data_io.SetFilePath(file_path);

	double A[2][3][2] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

	double_data_io.SetTotalNumber(12);

	double_data_io.WriteData((double*) A);

	double B[2][3][2];
	double_data_io.ReadData((double*) B);

	return 0;

}

#include <ctime>
#include <iostream>
#include <string>
#include <tuple>
#include<vector>
/*
 * Modified by Shiyong Tan, 2/7/18
 * Discard using matio, use DataIO instead
 * Start:
 */
//#include <matio.h>
#include "NumDataIO.h"
// End

//#include <C:/Program Files/MATLAB/R2016a/extern/include/mat.h> //Look for the path of mat.h
#include <Position.h>
#include <linterp.h>
#include <OTF.h>
#include "Common.h"

using namespace std;

// return an evenly spaced 1-corrSize grid of doubles.
std::vector<double> linspace(double first, double last, int len) {
	std::vector<double> result(len);
	double step = (last - first) / (len - 1);
	for (int i = 0; i<len; i++) { result[i] = first + i*step; }
	return result;
}

OTF::OTF(OTF& otf) {
	num_element = otf.GetNumElement();
	ncams = otf.GetNumCam();
	aData = new double*[ncams];
	bData = new double*[ncams];
	cData = new double*[ncams];
	alphaData = new double*[ncams];
	otf.GetAData(aData);
	otf.GetBData(bData);
	otf.GetCData(cData);
	otf.GetAlphaData(alphaData);
}

void OTF::GetAData(double** a_data) {
	for (int camid = 0; camid < ncams; camid++) {
		a_data[camid] = new double[num_element];
		for (int n = 0; n < num_element; n++) {
			a_data[camid][n] = aData[camid][n];
		}
	}
}

void OTF::GetBData(double** b_data) {
	for (int camid = 0; camid < ncams; camid++) {
		b_data[camid] = new double[num_element];
		for (int n = 0; n < num_element; n++) {
			b_data[camid][n] = bData[camid][n];
		}
	}
}

void OTF::GetCData(double** c_data) {
	for (int camid = 0; camid < ncams; camid++) {
		c_data[camid] = new double[num_element];
		for (int n = 0; n < num_element; n++) {
			c_data[camid][n] = cData[camid][n];
		}
	}
}

void OTF::GetAlphaData(double** alpha_data) {
	for (int camid = 0; camid < ncams; camid++) {
		alpha_data[camid] = new double[num_element];
		for (int n = 0; n < num_element; n++) {
			alpha_data[camid][n] = alphaData[camid][n];
		}
	}
}


//OTF::OTF(int cams, string& path) {
OTF::OTF(int ncams, string matfile) : ncams(ncams), mat_path(matfile) {

//	OTF_matdata();
	OTF_txtdata();
}

void OTF::OTF_txtdata() {
	// TODO Add error if no it does not find a text file
	// remove comments from the file 
	ifstream infile(mat_path.c_str(), ios::in);
	//ofstream outfile("parsed.txt");
	string line;
	aData = new double*[ncams];
	bData = new double*[ncams];
	cData = new double*[ncams];
	alphaData = new double*[ncams];

	stringstream parsed;

	while (getline(infile, line)) {
		size_t commentpos = line.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos) {
				line.erase(commentpos);
			}
			parsed << line << '\t';
		}
	}
	infile.close();

	int nElements; parsed >> nElements;
	num_element = nElements;
//	int nsubvector; parsed >> nsubvector;

	for (int camid = 0; camid < ncams; camid++) {
		aData[camid] = new double[nElements];
		bData[camid] = new double[nElements];
		cData[camid] = new double[nElements];
		alphaData[camid] = new double[nElements];
	}

	for (int camid = 0; camid < ncams; camid++) {

		for (int n = 0; n < nElements; n++) {
			parsed >> aData[camid][n];
		}
	}

	for (int camid = 0; camid < ncams; camid++) {

		for (int n = 0; n < nElements; n++) {
			parsed >> bData[camid][n];
		}
	}

	for (int camid = 0; camid < ncams; camid++) {

		for (int n = 0; n < nElements; n++) {
			parsed >> cData[camid][n];
		}
	}

	for (int camid = 0; camid < ncams; camid++) {

		for (int n = 0; n < nElements; n++) {
			parsed >> alphaData[camid][n];
		}
	}


//	for (int n = 0; n < nsubvector; n++) {
//		double x; parsed >> x;
//		gridx.push_back(x);
//	}
//
//	for (int n = 0; n < nsubvector; n++) {
//		double y; parsed >> y;
//		gridy.push_back(y);
//	}
//
//	for (int n = 0; n < nsubvector; n++) {
//		double z; parsed >> z;
//		gridz.push_back(z);
//	}

	/*outfile << parsed.rdbuf();
	outfile.close();*/
}
// obtaning the OTF parameter data from mat file
void OTF::OTF_matdata() {
//	const char *fileName = mat_path.c_str();
	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Discard using matio, using DataIO instead.
	 * Start:
	 */
//	mat_t *mat = Mat_Open(fileName, MAT_ACC_RDONLY);
//
//	if (mat == NULL) {
//		cout << " error in reading the mat file" << endl;
//	}
//	string cam;
//	string name;
//	string namesize;
//	string dataName;
//
//	aData = new double*[ncams];
//	bData = new double*[ncams];
//	cData = new double*[ncams];
//	alphaData = new double*[ncams];
//
//	for (int i = 0; i < ncams; i++) {
//
//		if (mat) {
//			//std::cout << "Open file to read\n\tmat == " << mat << "\n";
//
//			matvar_t *matVar = 0;
//			cam = to_string(i);
//			name = "A" + cam;
//			matVar = Mat_VarRead(mat, (char*)name.c_str());
//
//			if (matVar)
//			{
//				int rows;	int cols;
//				rows = matVar->dims[0]; cols = matVar->dims[1];
//				aData[i] = new double[rows*cols];
//				double *namedata = static_cast<double*>(matVar->data);
//				for (int j = 0; j < rows*cols; j++) {
//					aData[i][j] = namedata[j];
//				}
//			}
//
//			name = "B" + cam;
//			matVar = Mat_VarRead(mat, (char*)name.c_str());
//			if (matVar)
//			{
//				int rows;	int cols;
//				rows = matVar->dims[0]; cols = matVar->dims[1];
//				bData[i] = new double[rows*cols]; //(double *)malloc(sizeof(double)*rows*cols);
//				double *namedata = static_cast<double*>(matVar->data);
//				for (int j = 0; j < rows*cols; j++) {
//					bData[i][j] = namedata[j];
//				}
//			}
//
//			name = "C" + cam;
//			matVar = Mat_VarRead(mat, (char*)name.c_str());
//			if (matVar)
//			{
//				int rows;	int cols;
//				rows = matVar->dims[0]; cols = matVar->dims[1];
//
//				cData[i] = new double[rows*cols]; //(double *)malloc(sizeof(double)*rows*cols);
//				double *namedata = static_cast<double*>(matVar->data);
//				for (int j = 0; j < rows*cols; j++) {
//					cData[i][j] = namedata[j];
//				}
//			}
//
//			name = "Alpha" + cam;
//			matVar = Mat_VarRead(mat, (char*)name.c_str());
//			if (matVar)
//			{
//				int rows;	int cols;
//				rows = matVar->dims[0]; cols = matVar->dims[1];
//				namesize = "alpha" + cam + "Size";
//				unsigned namesize = matVar->nbytes / matVar->data_size;
//				dataName = "alpha" + cam + "Data";
//
//				alphaData[i] = new double[rows*cols]; //(double *)malloc(sizeof(double)*rows*cols);
//				double *namedata = static_cast<double*>(matVar->data);
//				for (int j = 0; j < rows*cols; j++) {
//					alphaData[i][j] = namedata[j];
//				}
//			}
//		}
//
//		else
//		{
//			cout << "Cannot open file\n";
//		}
//
//	}
//	Mat_Close(mat);
	// TODO: to check whether it works.
	aData = new double*[ncams];
	bData = new double*[ncams];
	cData = new double*[ncams];
	alphaData = new double*[ncams];

	// Data Format: A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4, Alpha1, Alpha2, Alpha3, Alpha4
	NumDataIO<double> data_io;
	data_io.SetFilePath(mat_path);
	int total_num = data_io.GetTotalNumber();
	double cam_data[total_num];
	data_io.SetTotalNumber(total_num);
	data_io.ReadData((double*) cam_data);

	for (int i = 0; i < ncams; i++) {
		int rows = total_num / 16;
		for (int j = 0; j < rows; j++ ) {
			aData[i][j] = cam_data[i * rows + j];
			bData[i][j] = cam_data[(i + 4) * rows + j];
			cData[i][j] = cam_data[(i + 8) * rows + j];
			alphaData[i][j] = cam_data[(i + 12) * rows + j];
		}
	}
	// End
}

//tuple<InterpMultilinear<3, double>, InterpMultilinear<3, double>, InterpMultilinear<3, double>, InterpMultilinear<3, double>> OTF::OTFgrid(int camid) {

vector <double> OTF::OTFgrid(int camID, Position& pos3D) {

	double cornerx1 = config.x_lower_limit, cornery1 = config.y_lower_limit, cornerz1 = config.z_lower_limit;
	double cornerx2 = config.x_upper_limt, cornery2 = config.y_upper_limt, cornerz2 = config.z_upper_limt;
	// construct the grid in each dimension. 
	// note that we will pass in a sequence of iterators pointing to the beginning of each grid

	array<int, 3> grid_sizes;
	//grid_sizes[0] = 3;
//	grid_sizes[1] = 3;
//	grid_sizes[2] = 3;

	grid_sizes[0] = (int)cbrt(num_element);
	grid_sizes[1] = grid_sizes[0];
	grid_sizes[2] = grid_sizes[0];


	std::vector<double> gridx = linspace(cornerx1, cornerx2, grid_sizes[0]);
	std::vector<double> gridy = linspace(cornery1, cornery2, grid_sizes[1]);
	std::vector<double> gridz = linspace(cornerz1, cornerz2, grid_sizes[2]);
	vector< vector<double>::iterator > grid_iter_list;
	grid_iter_list.push_back(gridx.begin());
	grid_iter_list.push_back(gridy.begin());
	grid_iter_list.push_back(gridz.begin());

	// the size of grid in each dimension
//	array<int, 3> grid_sizes;
//	//grid_sizes[0] = 3;
////	grid_sizes[1] = 3;
////	grid_sizes[2] = 3;
//
//	grid_sizes[0] = (int)pow(num_element, 1/3);
//	grid_sizes[1] = grid_sizes[0];
//	grid_sizes[2] = grid_sizes[0];

	// total number of elements
	int num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];

	// construct the interpolator, the last two arguments are pointers to the underlying data
//	double* test = aData[camID];
	InterpMultilinear<3, double> otfgrid_a(grid_iter_list.begin(), grid_sizes.begin(), aData[camID], aData[camID] + num_elements);
	InterpMultilinear<3, double> otfgrid_b(grid_iter_list.begin(), grid_sizes.begin(), bData[camID], bData[camID] + num_elements);
	InterpMultilinear<3, double> otfgrid_c(grid_iter_list.begin(), grid_sizes.begin(), cData[camID], cData[camID] + num_elements);
	InterpMultilinear<3, double> otfgrid_alpha(grid_iter_list.begin(), grid_sizes.begin(), alphaData[camID], alphaData[camID] + num_elements);

	vector<double> otfParam(4);
	array<double, 3> pos1 = { pos3D.X(), pos3D.Y(), pos3D.Z() };
	otfParam[0] = otfgrid_a.interp(pos1.begin());
	otfParam[1] = otfgrid_b.interp(pos1.begin());
	otfParam[2] = otfgrid_c.interp(pos1.begin());
	otfParam[3] = otfgrid_alpha.interp(pos1.begin());

	return otfParam;
}

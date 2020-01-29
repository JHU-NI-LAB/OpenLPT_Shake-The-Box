#ifndef OTF_H
#define	OTF_H

#include <ctime>
#include <iostream>
#include <string>
#include <assert.h>

//TODO: Temporary modification by Shiyong Tan, 1/31/18
//#include <matio.h>
#include <Position.h>

using namespace std;
class OTF {
public:

	//constructor: does nothing
	OTF() {};
	OTF(OTF& otf);
	//takes the number of cams and path to the .mat file (variable 'a' for cam0 should be saved as A0)
	//OTF(int cams, string& path);
	OTF(int cams, string matfile);
	//OTF(OTF& o);
	~OTF() {
		for (int n = 0; n < ncams; n++) {
			delete[] aData[n];
			delete[] bData[n];
			delete[] cData[n];
			delete[] alphaData[n];
		}
		delete[] aData;
		delete[] bData;
		delete[] cData;
		delete[] alphaData;
		//cout << "OTF destructor called" << endl;
	};

	// extract OTF values from .mat file
	void OTF_matdata();
	// extract OTF values from .txt file
	void OTF_txtdata();

	int GetNumElement() { return num_element; };
	int GetNumCam() { return ncams; };

	void GetAData(double** a_data);
	void GetBData(double** b_data);
	void GetCData(double** c_data);
	void GetAlphaData(double** alpha_data);

	// creating a grid for trilinear interpolation
	//tuple<InterpMultilinear<3, double>, InterpMultilinear<3, double>, InterpMultilinear<3, double>, InterpMultilinear<3, double>> OTF::OTFgrid(int camid);
	vector <double> OTFgrid(int camid, Position& pos3D);

	// gives OTF parameters at pos
	//vector <double> OTF::OTFParam(int camid, array<double, 3> pos);

private:
	int ncams;
	int num_element = 0;
	double **aData;
	double **bData;
	double **cData;
	double **alphaData;
	vector<double> gridx;
	vector<double> gridy;
	vector<double> gridz;
	string mat_path;
};

/*inline OTF::OTF(OTF& o) :
	ncams(o.ncams),
	aData(o.aData),
	bData(o.bData),
	cData(o.cData),
	alphaData(o.alphaData),
	mat_path(o.mat_path)
{}*/

#endif

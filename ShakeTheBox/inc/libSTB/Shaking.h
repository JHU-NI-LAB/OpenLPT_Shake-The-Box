#ifndef Shaking_H
#define	Shaking_H

#include <string>
#include <deque>

#include <Frame.h>
#include <Camera.h>
#include <OTF.h>

using namespace std;
class Shaking {
//class Shaking {
public:
	// ################################ STRUCTS #############################################
	// A structure to save pixel range for 1 x psize and 2 x psize
	struct PixelRange {
		PixelRange(int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4)
			: xmin1(x1), xmin2(x2), xmax1(x3), xmax2(x4), ymin1(y1), ymin2(y2), ymax1(y3), ymax2(y4)
		{}
		int xmin1, xmin2, xmax1, xmax2, ymin1, ymin2, ymax1, ymax2;
	};
	// ############################### FUNCTIONS ############################################

	// constructor
	Shaking(int ncams, int ignoreCam, OTF& otfcalib, int Npixw, int Npixh, int psize, double d, Position p, deque<Camera> camParam, deque<int**>& residual, double I);
	
	// destructor/home/shiyongtan/Documents/Research/AshwanthCode/Shake-the-box/CPVGDFTracker_gv/lib/OTF.cpp
	~Shaking() {
		for (int n = 0; n < rcams; n++) {
			for (int i = 0; i < 2 * psize; i++) {
//			for (int i = 0; i < psize; i++) {
				delete[] pixels_Part[n][i];
				delete[] pixels_PartAugRes[n][i];
			}
			delete[] pixels_Part[n];
			delete[] pixels_PartAugRes[n];
		}
		delete[] pixels_Part;
		delete[] pixels_PartAugRes;
		//cout << "Shaking destructor" << endl;
	};

	// gives the updated position after shaking
	void Pos();
	// gives the range of pixels on all cameras around the particle pos
	pair< deque<Shaking::PixelRange>, deque<Position> > PixRange(Position pos);
	// gives the intensity of reprojection at (x,y): (I_p(x,y))
	double PartReproj(Position particle2Dcenter, vector<double>& otfParam, int x, int y);
	// Particle Augemented Residual Matrix: (I_res+p)
	void PartAugResImage(int camID, int id, PixelRange p);
	// gives the residual of updated 3D position: R
	double Res(Position posnew);

	// creating a second order polynomial
	deque< double> Quadratic(deque<double> X, deque < double > R);

	// updates that particles 3D intensity
	void Int();

	// adds 2D position and intensity data to pos
	void FullData(Position& pos, double intensity, int Cams, int ignoreCam);

	// returns the index of smallest element of an array
	int IndexofSmallestElement(deque<double> array, int size);
	int IndexofLargestElement(double *array, int size);

	Position& Get_posnew();
	double Get_int();

//private:
	//################################## VARIABLES ##########################################

	int psize;
	
	Position pos3Dold;
	deque<Position> pos2Dold, pos2Dnew;
	deque<PixelRange> pRangeOld, pRangeNew;

	//original image
	deque<int**> pixels_orig;
	//residual image
	deque<int**> pixels_Res;
	//OTF Parameters
	OTF OTFcalib; vector<double> otfParamold;
	//reprojection
	int*** pixels_Part;
	//particle augmented reprojection
	int*** pixels_PartAugRes;
	//3D correction size
	double corrSize;
	//cam parameters
	deque<Camera> camParam;	int ncams;	int Npixh;	int Npixw; int rcams;
	int ignoreCam;
	// new position and intensity
	Position pos3Dnew;
	double int3D;
};

inline Position& Shaking::Get_posnew() {
	return Shaking::pos3Dnew;
}

inline double Shaking::Get_int() {
	return int3D;
}

#endif

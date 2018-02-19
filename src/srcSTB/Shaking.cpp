#include <string>
#include <deque>

#include <Shaking.h>

using namespace std;



// constructors
Shaking::Shaking(int Ncams, int ignoredCam, OTF& otfcalib, int NpixW, int NpixH, int Psize, double corrsize, Position p, deque<Camera> camparam, deque<int**>& residual, double I) :
	ncams(Ncams), ignoreCam(ignoredCam), Npixw(NpixW), Npixh(NpixH), psize(Psize), corrSize(corrsize), pos3Dold(p), camParam(camparam), OTFcalib(otfcalib), pixels_Res(residual), int3D(I) {

	// Particle augmented reprojection (for all cameras)
	if (ignoreCam < ncams) {	rcams = ncams-1;	}
	else {	rcams = ncams;	}

	pixels_PartAugRes = new int**[rcams];
	pixels_Part = new int**[rcams];

	for (int i = 0; i < rcams; i++) {
		pixels_PartAugRes[i] = new int*[2*psize];
		pixels_Part[i] = new int*[2 * psize];
		for (int j = 0; j < 2 * psize; j++) {
			pixels_PartAugRes[i][j] = new int[2*psize];
			pixels_Part[i][j] = new int[2 * psize];
		}
	}

	//Getting the pixel 2D centers and its range around the center on each camera
	pair< deque<PixelRange>, deque<Position> > pOld = PixRange(pos3Dold);
	pRangeOld = pOld.first; pos2Dold = pOld.second;

	// creating a particle reproj image (I_p) matrix in the pixel range
	int id = 0;
	for (int camID = 0; camID < ncams; camID++) {
		if (camID != ignoreCam) {
			// OTF parameters at old 3D position
			otfParamold = OTFcalib.OTFgrid(camID, pos3Dold);
			for (int x = pRangeOld[id].xmin2; x < pRangeOld[id].xmax2; x++) {
				for (int y = pRangeOld[id].ymin2; y < pRangeOld[id].ymax2; y++) {
					int X = x - pRangeOld[id].xmin2; int Y = y - pRangeOld[id].ymin2;
					pixels_Part[id][Y][X] = round(PartReproj(pos2Dold[id], otfParamold, x, y));
				}
			}
			id++;
		}
	}

	// creating a particle augmented residual image: (I_(res+p))
	id = 0;
	for (int camID = 0; camID < ncams; camID++) {
		if (camID != ignoreCam) {
			PartAugResImage(camID, id, pRangeOld[id]);
			id++;
		}
	}
	
	//updating the particle position and intensity
	Pos();
	// pixel range on the cameras around the updated particle centers
	pair< deque<PixelRange>, deque<Position> > pNew = PixRange(pos3Dnew);
	pRangeNew = pNew.first; pos2Dnew = pNew.second;
	// updating the 3D intensity
	Int();
}

// updates the 3D position by shaking
void Shaking::Pos() {

	/*
	 * Modified by Shiyong Tan, 1/31/18
	 * Illegal initialization with {...}
	 * Start:
	 */
//	deque<double> del = { -corrSize, 0, +corrSize }; // mm correction in 3D
	deque<double> del(3);
	del[0] = -corrSize; del[1] = 0; del[2] = +corrSize;
	// END
	pos3Dnew = pos3Dold;
	deque < double > R(4); //saves the residual at -del, 0, +del and local minima
	int size;
															// for X first
	deque<double> X(4), Y(4), Z(4);
	for (int i = 0; i < 3; i++) {
		X[i] = (pos3Dold).X() + del[i];
		pos3Dnew.Set_X(X[i]); // shifting by del in X-direction
		R[i] = Res(pos3Dnew);
	}
	
	deque<double> coeffx = Quadratic(X, R);	// finding the coefficients of a quadratic that passes through the three residual values
	X[3] = -coeffx[1] / (2 * coeffx[0]); // finding the minima of the quadratic 

	// consider the minima only if it lies between x-del and x+del
	if (X[3] > X[0] && X[3] < X[2]) {
		size = 4; 
		pos3Dnew.Set_X(X[3]); R[3] = Res(pos3Dnew); // calulating the actual residual at minima
	}
	else
		size = 3;

	int x3D = IndexofSmallestElement(R, size); pos3Dnew.Set_X(X[x3D]);	// setting X to the position with smallest residual

															// then update Y
	for (int i = 0; i < 3; i++) {
		Y[i] = (pos3Dold).Y() + del[i];
		pos3Dnew.Set_Y(Y[i]);
		R[i] = Res(pos3Dnew);
	}
	deque<double> coeffy = Quadratic(Y, R);
	Y[3] = -coeffy[1] / (2 * coeffy[0]); 

	if (Y[3] > Y[0] && Y[3] < Y[2]) {
		size = 4; 
		pos3Dnew.Set_Y(Y[3]); R[3] = Res(pos3Dnew); // calulating the actual residual at minima
	}
	else
		size = 3;
	int y3D = IndexofSmallestElement(R, size); pos3Dnew.Set_Y(Y[y3D]);
	
															// then update Z
	for (int i = 0; i < 3; i++) {
		Z[i] = (pos3Dold).Z() + del[i];
		pos3Dnew.Set_Z(Z[i]);
		R[i] = Res(pos3Dnew);
	}
	deque<double> coeffz = Quadratic(Z, R);
	Z[3] = -coeffz[1] / (2 * coeffz[0]); 

	if (Z[3] > Z[0] && Z[3] < Z[2]) {
		size = 4; 
		pos3Dnew.Set_Z(Z[3]); R[3] = Res(pos3Dnew); // calulating the actual residual at minima
	}
	else
		size = 3;
	int z3D = IndexofSmallestElement(R, size); pos3Dnew.Set_Z(Z[z3D]);
}

// takes the new 3D particle position and gives the residual R
double Shaking::Res(Position posnew) {
	deque<Position> particle2Dnew;
	double R = 0;
	int id = 0;
	for (int camID = 0; camID < ncams; camID++) {
		if (camID != ignoreCam) {
			//old and new 2D projection centers
			Position pnew(camParam[camID].WorldToImage(posnew)); pnew.Set_Z(0);
			particle2Dnew.push_back(camParam[camID].Distort(pnew));

			//otf parameters at new position
			vector <double> otfParamnew = OTFcalib.OTFgrid(camID, posnew);

			// Calculating the residual for updated 3D position
			for (int x = pRangeOld[id].xmin2; x < pRangeOld[id].xmax2; x++) {
				for (int y = pRangeOld[id].ymin2; y < pRangeOld[id].ymax2; y++) {
					int X = x - pRangeOld[id].xmin2, Y = y - pRangeOld[id].ymin2;
					R = R + pow((pixels_PartAugRes[id][Y][X] - round(PartReproj(particle2Dnew[id], otfParamnew, x, y))), 2);
				}
			}
			id++;
		}
	}
	return R;
}

// gives the range of pixels at psize and 2*psize around the 2D centers (on each camera)
pair< deque<Shaking::PixelRange>, deque<Position> > Shaking::PixRange(Position pos) {
	deque<PixelRange> pRange;
	deque<Position> pos2D;
	int id = 0;
	for (int camID = 0; camID < ncams; camID++) {
		if (camID != ignoreCam) {
			Position p(camParam[camID].WorldToImage(pos)); p.Set_Z(0);
			pos2D.push_back(camParam[camID].Distort(p));
			// pixel range for each particle with 1 x particle size and 2 x particle size
			int xmin1 = max(1, (int)floor(pos2D[id].X() - psize / 2));	 int xmin2 = max(1, (int)floor(xmin1 - psize / 2));
			int ymin1 = max(1, (int)floor(pos2D[id].Y() - psize / 2));	 int ymin2 = max(1, (int)floor(ymin1 - psize / 2));
			int xmax1 = min(Npixw, (int)floor(pos2D[id].X() + psize / 2)); int xmax2 = min(Npixw, (int)ceil(xmax1 + psize / 2));
			int ymax1 = min(Npixh, (int)floor(pos2D[id].Y() + psize / 2)); int ymax2 = min(Npixh, (int)ceil(ymax1 + psize / 2));
			PixelRange pixRange(xmin1, xmin2, xmax1, xmax2, ymin1, ymin2, ymax1, ymax2);
			pRange.push_back(pixRange);
			id++;
		}
	}

	return make_pair(pRange, pos2D);
}

// a function for Gaussian ellipse reprojection at position (x,y)
double Shaking::PartReproj(Position particle2Dcenter, vector <double>& otfParam, int x, int y) {
	double xx = ((double)x - particle2Dcenter.X())*cos(otfParam[3]) + ((double)y - particle2Dcenter.Y())*sin(otfParam[3]);
	double yy = -((double)x - particle2Dcenter.X())*sin(otfParam[3]) + ((double)y - particle2Dcenter.Y())*cos(otfParam[3]);
	double value = otfParam[0] * exp(-(otfParam[1] * pow(xx, 2) + otfParam[2] * pow(yy, 2)));
	return(value);
}

// takes the camID, 2D projection center and gives particle augmented residual on that camera
void Shaking::PartAugResImage(int camID, int id, PixelRange p) {
	for (int x = p.xmin2; x < p.xmax2; x++) {
		for (int y = p.ymin2; y < p.ymax2; y++) {
			int X = x - p.xmin2, Y = y - p.ymin2;
			pixels_PartAugRes[id][Y][X] = max(0,min(255,pixels_Res[camID][y][x] + pixels_Part[id][Y][X]));
		}		      //aug intensity   =  //residual intensity   +    //particle intensity
	}
}

// returns the coefficients of quadratic passing through 3 points
deque<double> Shaking::Quadratic(deque<double> X, deque<double> R) {
	// R = a*x^2 + b*x + c
	deque <double> coeff(3);
	coeff[0] = (((R[1] - R[0]) / (X[1] - X[0])) - ((R[2] - R[0]) / (X[2] - X[0]))) / (X[1] - X[2]); //a
	coeff[1] = ((R[1] - R[0]) / (X[1] - X[0])) - coeff[0]*(X[1] + X[0]); //b
	coeff[2] = R[2] - coeff[0] * pow(X[2], 2) - coeff[1]*X[2]; //c

	return coeff;

}

// updates the intesity of particle
void Shaking::Int() {
	double num = 0.0, denum = 0.0;
	int id = 0;
	double* peakIntensity = new double[rcams];
	for (int i = 0; i < rcams; i++)
		peakIntensity[i] = 0.0;
	deque<Position> pos2D;
	deque< vector<double> > otfparam;

	// finding the camera with highest local peak intensity
	for (int camID = 0; camID < ncams; camID++) {
		if (camID != ignoreCam) {
			otfparam.push_back(OTFcalib.OTFgrid(camID, pos3Dnew));
			int xmin = max(pRangeOld[id].xmin1, pRangeNew[id].xmin1), xmax = min(pRangeNew[id].xmax1, pRangeOld[id].xmax1);
			int ymin = max(pRangeOld[id].ymin1, pRangeNew[id].ymin1), ymax = min(pRangeNew[id].ymax1, pRangeOld[id].ymax1);
			for (int x = xmin; x < xmax; x++) {
				for (int y = ymin; y < ymax; y++) {
					int X = x - pRangeOld[id].xmin2, Y = y - pRangeOld[id].ymin2;
					peakIntensity[id] = max(peakIntensity[id], PartReproj(pos2Dnew[id], otfparam[id], x, y));
				}
			}
			id++;
		}
	}

	// updating particle intensity (ignoring the camera with highest local peak intensity)
	double ignore = IndexofLargestElement(peakIntensity, rcams);
	for (int ID = 0; ID < rcams; ID++) {
		if (ID != ignore) {
			int xmin = max(pRangeOld[ID].xmin1, pRangeNew[ID].xmin1), xmax = min(pRangeNew[ID].xmax1, pRangeOld[ID].xmax1);
			int ymin = max(pRangeOld[ID].ymin1, pRangeNew[ID].ymin1), ymax = min(pRangeNew[ID].ymax1, pRangeOld[ID].ymax1);
			for (int x = xmin; x < xmax; x++) {
				for (int y = ymin; y < ymax; y++) {
					int X = x - pRangeOld[ID].xmin2, Y = y - pRangeOld[ID].ymin2;
					num = num + pixels_PartAugRes[ID][Y][X];
					denum = denum + PartReproj(pos2Dnew[ID], otfparam[ID], x, y);
				}
			}
		}
	}
	delete[] peakIntensity;
	//cout << ignore << "," << num << "," << denum << endl;
	int3D = int3D * sqrt(abs(num / denum));
}

// returns index of smallest element
int Shaking::IndexofSmallestElement(deque<double> array, int size) {
	int index = 0;
	for (int i = 1; i < size; i++) {
		if (array[i] < array[index])
			index = i;
	}
	return index;
}

// returns index of largest element
int Shaking::IndexofLargestElement(double *array, int size) {
	int index = 0;
	for (int i = 1; i < size; i++) {
		if (array[i] > array[index])
			index = i;
	}
	return index;
}


#include <string>
#include <deque>

#include <Shaking.h>
#include "NumDataIO.h"

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
//	for (int i = 0; i < rcams; i++) {
//		pixels_PartAugRes[i] = new int*[psize];
//		pixels_Part[i] = new int*[psize];
//		for (int j = 0; j < psize; j++) {
//			pixels_PartAugRes[i][j] = new int[psize];
//			pixels_Part[i][j] = new int[psize];
//			for (int k = 0; k < psize; k++) {
//				pixels_Part[i][j][k] = 0;
//				pixels_PartAugRes[i][j][k] = 0; //initialize
//			}
//		}
//	}

	//Getting the pixel 2D centers and its range around the center on each camera
	pair< deque<PixelRange>, deque<Position> > pOld = PixRange(pos3Dold);
	pRangeOld = pOld.first; pos2Dold = pOld.second;

	// creating a particle reproj image (I_p) matrix in the pixel range
	int id = 0;
	for (int camID = 0; camID < ncams; camID++) {
		if (camID != ignoreCam) {
			// OTF parameters at old 3D position
			otfParamold = OTFcalib.OTFgrid(camID, pos3Dold);
//			for (int x = pRangeOld[id].xmin1; x < pRangeOld[id].xmax1; x++) {
//				for (int y = pRangeOld[id].ymin1; y < pRangeOld[id].ymax1; y++) {
//					int X = x - pRangeOld[id].xmin1; int Y = y - pRangeOld[id].ymin1;
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
	
	// output the residual image
//	NumDataIO<int> image_output;
//	string save_path;
//	for (int n = 0; n < ncams; n++) {
//		save_path = "/home/sut210/Documents/Experiment/EXP6/Augresimg" + to_string(n) +".txt";
//		image_output.SetFilePath(save_path);
//		image_output.SetTotalNumber(Npixh * Npixw);
//		int* pixel_array = new int[Npixh * Npixw];
//		for (int i = 0; i < Npixh; i++)
//			for (int j = 0; j < Npixw; j++)
//				pixel_array[i * Npixw + j] = pixels_PartAugRes[n][i][j];
//		image_output.WriteData(pixel_array);
//	}

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

			//NumDataIO<int> data_io;
			//int* pixels_partaugres;
			//int* part_proj;
			//if (pRangeOld[id].xmax2 > pRangeOld[id].xmin2 && pRangeOld[id].ymax2 > pRangeOld[id].ymin2) {
			//	int num = (pRangeOld[id].ymax2 - pRangeOld[id].ymin2) * (pRangeOld[id].xmax2 - pRangeOld[id].xmin2 );
			//	pixels_partaugres = new int[num];
			//	part_proj= new int[num];
			//}
			pair< deque<PixelRange>, deque<Position> > pNew = PixRange(posnew);
			deque<PixelRange> pNew_range = pNew.first;

			// Calculating the residual for updated 3D position
			for (int x = pRangeOld[id].xmin2; x < pRangeOld[id].xmax2; x++) {
				for (int y = pRangeOld[id].ymin2; y < pRangeOld[id].ymax2; y++) {
//			for (int x = pRangeOld[id].xmin1; x < pRangeOld[id].xmax1; x++) {
//				for (int y = pRangeOld[id].ymin1; y < pRangeOld[id].ymax1; y++) {
//			for (int x = pNew_range[id].xmin1; x < pNew_range[id].xmax1; x++) {
//				for (int y = pNew_range[id].ymin1; y < pNew_range[id].ymax1; y++) {
					int X = x - pRangeOld[id].xmin2, Y = y - pRangeOld[id].ymin2;
					if (Y < 2 * psize && Y >= 0 && X < 2 * psize && X >= 0) {
						R = R + pow((pixels_PartAugRes[id][Y][X] - round(PartReproj(particle2Dnew[id], otfParamnew, x, y))), 2);
					} else {
						R = R + pow((pixels_Res[id][y][x] - round(PartReproj(particle2Dnew[id], otfParamnew, x, y))), 2);
					}
					//pixels_partaugres[X * (pRangeOld[id].ymax2 - pRangeOld[id].ymin2) + Y] = pixels_PartAugRes[id][Y][X];
					//part_proj[X * (pRangeOld[id].ymax2 - pRangeOld[id].ymin2) + Y] = round(PartReproj(particle2Dnew[id], otfParamnew, x, y));
				}
			}
			//if (pRangeOld[id].xmax2 > pRangeOld[id].xmin2 && pRangeOld[id].ymax2 > pRangeOld[id].ymin2) {
			//	int num = (pRangeOld[id].ymax2 - pRangeOld[id].ymin2) * (pRangeOld[id].xmax2 - pRangeOld[id].xmin2 );
			//data_io.SetTotalNumber(num);
			//data_io.SetFilePath("/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/Projection/pixels_partaugres.txt");
			//data_io.WriteData((int*) pixels_partaugres);
			//data_io.SetFilePath("/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/Projection/part_proj.txt");
			//data_io.WriteData((int*) part_proj);
			//delete[] pixels_partaugres;
			//delete[] part_proj;
			//}
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
//			int xmin1 = max(1, (int)floor(pos2D[id].X() - psize / 4));	 int xmin2 = max(1, (int)floor(xmin1 - psize / 4));
//			int ymin1 = max(1, (int)floor(pos2D[id].Y() - psize / 4));	 int ymin2 = max(1, (int)floor(ymin1 - psize / 4));
//			int xmax1 = min(Npixw, (int)floor(pos2D[id].X() + psize / 4)); int xmax2 = min(Npixw, (int)ceil(xmax1 + psize / 4));
//			int ymax1 = min(Npixh, (int)floor(pos2D[id].Y() + psize / 4)); int ymax2 = min(Npixh, (int)ceil(ymax1 + psize / 4));
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
	//NumDataIO<int> data_io;
	//int num = (p.xmax2 - p.xmin2) * (p.ymax2 - p.ymin2);
	//int* particle_res;
	//int* part_proj;
	//if (num > 0) {
	//	particle_res = new int[num];
	//	part_proj = new int[num];
	//}
	for (int x = p.xmin2; x < p.xmax2; x++) {
		for (int y = p.ymin2; y < p.ymax2; y++) {
	int X = x - p.xmin2, Y = y - p.ymin2;
//	for (int x = p.xmin1; x < p.xmax1; x++) {
//		for (int y = p.ymin1; y < p.ymax1; y++) {
//			int X = x - p.xmin1, Y = y - p.ymin1;
//			cout<<x<<","<<y<<",";
//			cout<<pixels_Res[camID][y][x]<<",";
//			cout<<pixels_Part[id][Y][X]<<",";
			pixels_PartAugRes[id][Y][X] = max(0,min(255,pixels_Res[camID][y][x] + pixels_Part[id][Y][X]));
//			particle_res[X * (p.ymax2 - p.ymin2) + Y] = pixels_Res[camID][y][x];
//			part_proj[X * (p.ymax2 - p.ymin2) + Y] = pixels_Part[id][Y][X];
//			cout<<pixels_PartAugRes[i[X * (p.ymax2 - p.ymin2) + Y]d][Y][X]<<endl;
		}		      //aug intensity   =  //residual intensity   +    //particle intensity
	}
//	if (num > 0) {
//		data_io.SetTotalNumber(num);
//		data_io.SetFilePath("/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/Projection/particle_res.txt");
//		data_io.WriteData(particle_res);
//		data_io.SetFilePath("/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/Projection/particle_proj.txt");
//		data_io.WriteData(part_proj);
//		delete[] particle_res;
//		delete[] part_proj;
//	}

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
	int* cam_ID = new int[rcams];
	double* peakIntensity = new double[rcams];
	for (int i = 0; i < rcams; i++)
		peakIntensity[i] = 0.0;
	deque<Position> pos2D;
	deque< vector<double> > otfparam;

	// finding the camera with highest local peak intensity
	for (int camID = 0; camID < ncams; camID++) {
		if (camID != ignoreCam) {
			cam_ID[id] = camID;
			otfparam.push_back(OTFcalib.OTFgrid(camID, pos3Dnew));
			int xmin = max(pRangeOld[id].xmin1, pRangeNew[id].xmin1), xmax = min(pRangeNew[id].xmax1, pRangeOld[id].xmax1);
			int ymin = max(pRangeOld[id].ymin1, pRangeNew[id].ymin1), ymax = min(pRangeNew[id].ymax1, pRangeOld[id].ymax1);
			for (int x = xmin; x < xmax; x++) {
				for (int y = ymin; y < ymax; y++) {
					int X = x - pRangeOld[id].xmin2, Y = y - pRangeOld[id].ymin2;

//					int X = x - pRangeOld[id].xmin1, Y = y - pRangeOld[id].ymin1;
//					peakIntensity[id] = max(peakIntensity[id], PartReproj(pos2Dnew[id], otfparam[id], x, y));
//					if (x < 1 || x > Npixw || y < 1 || y > Npixh) continue;
					if (Y < 2 * psize && Y >= 0 && X < 2 * psize && X >= 0) {
						peakIntensity[id] = max(peakIntensity[id], (double)pixels_PartAugRes[id][Y][X]);
					}
				}
			}
			id++;
		}
	}

	// updating particle intensity (ignoring the camera with highest local peak intensity)
	double ignore = IndexofLargestElement(peakIntensity, rcams);
//	double ratio = 255;
//	int nonzero_proj = 0;
	for (int ID = 0; ID < rcams; ID++) {
//		int num_single = 0;
		if (ID != ignore) {
//			int xmin = max(pRangeOld[ID].xmin1, pRangeNew[ID].xmin1), xmax = min(pRangeNew[ID].xmax1, pRangeOld[ID].xmax1);
//			int ymin = max(pRangeOld[ID].ymin1, pRangeNew[ID].ymin1), ymax = min(pRangeNew[ID].ymax1, pRangeOld[ID].ymax1);
			int xmin = pRangeNew[ID].xmin1, xmax = pRangeNew[ID].xmax1;
			int ymin = pRangeNew[ID].ymin1, ymax = pRangeNew[ID].ymax1;

			for (int x = xmin; x < xmax; x++) {
				for (int y = ymin; y < ymax; y++) {
					int X = x - pRangeOld[ID].xmin2, Y = y - pRangeOld[ID].ymin2;
//					int X = x - pRangeOld[ID].xmin1, Y = y - pRangeOld[ID].ymin1;
//					int X = x - xmin, Y = y - ymin;
					if (Y < 2 * psize && Y >= 0 && X < 2 * psize && X >= 0) { // the range of calculated PartAugRes
//					if (Y < psize && Y >= 0 && X < psize && X >= 0) {
						num = num + pixels_PartAugRes[ID][Y][X];
//						num_single = num_single + pixels_PartAugRes[ID][Y][X];
					} else {
						num = num + pixels_Res[cam_ID[ID]][y][x];
//						num_single = num_single + pixels_PartAugRes[ID][Y][X];
					}
					int pixel_proj = round(PartReproj(pos2Dnew[ID], otfparam[ID], x, y));
					denum = denum + pixel_proj;
//					if (pixel_proj > 0 && pixels_PartAugRes[ID][Y][X] / pixel_proj < ratio) {
//						ratio = pixels_PartAugRes[ID][Y][X] / pixel_proj;
//					}
				}
			}
//			if (num_single > 30) nonzero_proj ++;
//			int centerx = round((pRangeNew[ID].xmin1 + pRangeNew[ID].xmax1) / 2);
//			int centery = round((pRangeNew[ID].ymin1 + pRangeNew[ID].ymax1) / 2);
//			int proj_center = round(PartReproj(pos2Dnew[ID], otfparam[ID], centerx, centery));
//			if (proj_center > 0) {
//				int X = centerx - pRangeOld[ID].xmin2, Y = centery - pRangeOld[ID].ymin2;
//				double r = 255;
//				if (Y < psize * 2 && Y >= 0 && X < psize * 2 && X >= 0) {
//					r = pixels_PartAugRes[ID][Y][X] / proj_center;
//				} else {
//					r = pixels_Res[ID][centery][centerx] / proj_center;
//				}
//				if (r < ratio) {
//					ratio = r;
//				}
//			}
		}
	}


	// For reduced camera, if the particles are within the image, then add the intensity
	if (ignoreCam < ncams) {
		Position p(camParam[ignoreCam].WorldToImage(pos3Dnew)); p.Set_Z(0);
		Position pos2D_reducedcam = camParam[ignoreCam].Distort(p);
		// pixel range for each particle with 1 x particle size
		int xmin1 = (int)floor(pos2D_reducedcam.X() - psize / 2);
		//int xmin1 = max(1, (int)floor(pos2D_reducedcam.X() - psize / 2));
		int ymin1 = (int)floor(pos2D_reducedcam.Y() - psize / 2);
		//int ymin1 = max(1, (int)floor(pos2D_reducedcam.Y() - psize / 2));
		int xmax1 = (int)floor(pos2D_reducedcam.X() + psize / 2);
//		int xmax1 = min(Npixw, (int)floor(pos2D_reducedcam.X() + psize / 2));
		int ymax1 = (int)floor(pos2D_reducedcam.Y() + psize / 2);
//		int ymax1 = min(Npixh, (int)floor(pos2D_reducedcam.Y() + psize / 2));

//		if (xmin1 < Npixw && xmax1 > 1 && ymin1 < Npixh && ymax1 > 1 &&
//				pos2D_reducedcam.X() > 1 && pos2D_reducedcam.X() < Npixw &&
//				pos2D_reducedcam.Y() > 1 && pos2D_reducedcam.Y() < Npixh) {
		if (xmin1 > 1 && xmin1 < Npixw - psize  && xmax1 > psize && xmax1 < Npixw &&
				ymin1 > 1 && ymin1 < Npixh - psize  && ymax1 > psize && ymax1 < Npixh) {
			for (int x = xmin1; x < xmax1; x++) {
				for (int y = ymin1; y < ymax1; y++) {
					num = num + pixels_Res[ignoreCam][y][x];  // the residual for ignorecam is stored at the last.
					vector<double> OTFparam = OTFcalib.OTFgrid(ignoreCam, pos3Dnew);
					int pixel_proj = round(PartReproj(pos2D_reducedcam, OTFparam, x, y));
					denum = denum + pixel_proj;
				}
			}
		}

	}

	delete[] peakIntensity;
	delete[] cam_ID;
	//cout << ignore << "," << num << "," << denum << endl;
//	if (int3D == 0) int3D = 1;
	int3D = int3D * sqrt(abs(num / denum));// * ratio;
//	if (nonzero_proj < 2) {  // change to 2 because one camera with highest intensity has been ignore.
//		int3D = 0;  // indicate this particle disappears in at least two cameras.
//	}
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


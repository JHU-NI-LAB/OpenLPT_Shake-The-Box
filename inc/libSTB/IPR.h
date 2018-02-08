#ifndef IPR_H
#define	IPR_H

#include <string>
#include <deque>
#include <stdexcept>
#include <utility>

/* Modified by Shiyong Tan, 1/31/18
 * Supposed that it used the wrong head name, I think it is not needed.
 * Delete it directly. TODO: Check whether it counts.
 */
//#include <Tiff.h>
#include <Tiff2DFinder.h>
#include <Calibration.h>
#include <Shaking.h>

class IPR {

public:

	friend class STB;
	// default constructor: do nothing
	IPR() {};
	//constructor: file contains information about image names and OTF param
	IPR(std::string& fname, int ncams);
	//destructor
	~IPR() {
		for (int n = 0; n < ncams; n++) {
			if (pixels_orig.size() > 0) {
				for (int i = 0; i < Npixh; i++) {
					delete[] pixels_orig[n][i];
					delete[] pixels_reproj[n][i];
					delete[] pixels_res[n][i];
				}
				delete[] pixels_orig[n];
				delete[] pixels_reproj[n];
				delete[] pixels_res[n];
			}
		}
	};

	Frame FindPos3D(deque< deque<string> > filename, int frameNumber);
	Frame IPRLoop(Calibration& calib, OTF& OTFcalib, deque<int> camNums, int ignoreCam, double colors,
						deque<int**>& orig, deque<int**>& reproj, deque<int**>& res);

	// a fuction that takes all the 3D positions and gives reprojected images: (I_proj)
	void ReprojImage(Frame matched3D, OTF& OTFcalib, deque<int**>& pixels_reproj, bool STB);

	// gives the intensity of reprojection at pixel (x,y): (I_p(x,y))
	double PixelReproj(Position& particle2Dcenter, vector<double>& otfParam, int x, int y);

	// remove ghost and ambiguous paticles
	pair<int, int> Rem(Frame& pos3D, deque<double>& int3D, double mindist_3D);

	// adding the intesity and 2D position data
	void FullData(Position& pos, double intensity, int cams, int ignoreCam);

	// creating .mat files
	void MatfilePositions(deque<Position> pos3D, string mat_name);
	void MatfileImage(deque<int**>& pix, string name);
	void Load_2Dpoints(string name, int frame, int ignoreCam);

	// Get functions							
	deque<Position> Get_IPRpos();		//Get the 3D position and its intensity
	int Get_threshold();				// return the pixel threshold to find 2D particles
	double Get_mindist2D();				
	double Get_mindist3D();
	double Get_psize();
	double Get_intensityLower();
	string& Get_calibfile();
	bool triangulationOnly;
	bool IPROnly;

protected:
	int frame;
	deque<string> filename;

	// addresses
	string tiffaddress;
	string calibfile;
	string otfFile;
		
	// camera parameters
	std::deque<Camera> camsAll, cams4;
	std::deque<Camera> camsReduced;
	int ncams, Npixh, Npixw;
	bool reducedCams;

	// customizable input parameteres
	int it_outerloop, it_innerloop, it_reducedCam;
	int threshold, pixelMax;
	int psize;
	
	// OTF parameters
	//OTF OTFcalib;

	// 2D positions
	deque<Frame> iframes;

	//original images
	deque<int**> pixels_orig;  
	//repojected images
	deque<int**> pixels_reproj;  
	//residual image
	deque<int**> pixels_res;

	//tolerences
	double mindist_2D, mindist_3D;
	double mindist_2Dr, mindist_3Dr;

	//for ghost particle
	double intensityLower;

	// output: particle 3D positions and intensity
	deque<Position> pos3D;
	deque<double> intensity3D;

	// for testing:
	deque<Position> ghost3D;

private:
	void Position2Array(deque<Position> pos, double array[][12]);
};


deque<Position> inline IPR::Get_IPRpos() {
	return pos3D;
}

int inline IPR::Get_threshold() {
	return threshold;
}

double inline IPR::Get_mindist2D() {
	return mindist_2D;
}

double inline IPR::Get_mindist3D() {
	return mindist_3D;
}

double inline IPR::Get_psize() {
	return psize;
}

double inline IPR::Get_intensityLower() {
	return intensityLower;
}

inline string& IPR::Get_calibfile() {
	return calibfile;
}
#endif

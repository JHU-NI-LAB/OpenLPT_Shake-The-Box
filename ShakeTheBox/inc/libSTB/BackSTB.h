#ifndef BackSTB_H
#define	BackSTB_H

#include "STB.h"

using namespace std;

/*
 * Modified by Shiyong Tan, 2/1/18
 * STB.h includes this file, but this file also includes STB.
 * So, we have to declare STB here since this file is compiled ahead of STB.
 * Start:
 */
//class STB;
// End

class BackSTB {
public:

	// constructor
	BackSTB(int firstFrame, int lastFrame, double threshold);
	// destructor
	~BackSTB() {
		for (int n = 0; n < ncams; n++) {
			if (pixels_orig[n] != NULL)
				for (int i = 0; i < Npixh; i++) {
					delete[] pixels_orig[n][i];
					delete[] pixels_reproj[n][i];
					delete[] pixels_res[n][i];
				}
			delete[] pixels_orig[n];
			delete[] pixels_reproj[n];
			delete[] pixels_res[n];
		}
	};

	//############################### FUNCTIONS ##############################
	void UpdateTracks(STB& s);
	int Predictor(STB& s, int prevFrame, deque<Track>::iterator& Ltr,
		deque<Track>& unlinkedTracks, Frame& predictions,
		deque<double>& intensity);
	double LMSWienerPredictor(Track tracks, string direction, int order);
	vector<double> Polyfit(Track tracks, string direction, int prevFrame);
	void Shake(STB& s, Frame& estimate, deque<double>& intensity, deque<Track>& unlinkedTracks);
	deque<int> Rem(STB& s,Frame& pos3D, deque<double>& int3D, double mindist_3D, deque<Track>& unlinkedTracks);
	void MakeShortLinkResidual_backward(STB& s, int prevFrame, Frame& candidates, deque<Track>::iterator& tr, int iterations);
	
private:
	double distThreshold;										// threshold for connecting the track prediction to actual particle on a track
	int first, last;

	// camera parameters
	int ncams, Npixh, Npixw;

	// images
	deque<int**> pixels_orig;									// original images
	deque<int**> pixels_reproj;									// repojected images
	deque<int**> pixels_res;									// residual image
};


#endif

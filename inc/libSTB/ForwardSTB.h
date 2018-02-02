#ifndef ForwardSTB_H
#define	ForwardSTB_H

#include <STB.h>

using namespace std;

class ForwardSTB {
public:

	// constructor
	ForwardSTB(int firstFrame, int lastFrame, double threshold);
	// destructor
	~ForwardSTB() {
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
	void Predictor(STB& s, int prevFrame, deque<Track>::iterator& Ltr,
		deque<deque<Track>::iterator>& unlinkedTracks, Frame& predictions,
		deque<double>& intensity);
	void Shake(STB& s, Frame& estimate, deque<double>& intensity, deque<deque<Track>::iterator>& unlinkedTracks);
	deque<int> Rem(STB& s, Frame& pos3D, deque<double>& int3D, double mindist_3D, deque<deque<Track>::iterator>& unlinkedTracks);

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

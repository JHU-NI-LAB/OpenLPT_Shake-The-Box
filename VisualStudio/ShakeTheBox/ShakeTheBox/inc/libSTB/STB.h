#ifndef STB_H
#define	STB_H


#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <deque>
#include <stdexcept>
#include <utility>
#include <sstream>
#include <iomanip>

#include <Track.h>
#include <IPR.h>
#include <PredictiveField.h>
//#include <BackSTB.h>
//#include <ForwardSTB.h>

using namespace std;

class STB {
public:

//	friend class BackSTB;
//	friend class ForwardSTB;
	// default constructor
	STB() {};
	// constructor
	STB(int first, int last, string pfieldfile, string iprfile, int ncams, deque<int> camIDs, deque<string> imgNameFiles,
		double InitialRadius, double avgDist, double expShift, double masc, double mrsc, double fpt, double lowerInt, bool iprFlag);
	// copy constructor
	STB(STB& s);
	// destructor
	~STB() {
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

    enum TrackType{ Inactive = 0, ActiveShort = 1, ActiveLong = 2, Exit = 3, InactiveLong = 4, Buffer = 5};

	//############################### FUNCTIONS ##############################
	// to make tracks for the first four frames
	void InitialPhase(string pfieldfile);

	// a function to start a track for particles that were left untracked in current frame
	void StartTrack(int frame, PredictiveField& pField);
	// extends the particle track to nextFrame using search radius method (can be used for first 4 links in both, intialization and convergence phase)
	void MakeLink(int nextFrame, const Position& new_velocity, double radius, Track& track, bool& activeTrack);
	pair<Frame::const_iterator, float> ComputeCost(Frame& fr1, Frame& fr2, double radius,
		const Position& estimate,
		const Position& velocity,
		const Position& new_velocity,
		bool stopflag);
	pair<Frame::const_iterator, float> ComputeCost(Frame& fr1, Frame& fr2, double radius,
		const Position& estimate,
		const Position& velocity,
		const Position& new_velocity,
		bool stopflag,
		bool* candidate_used);

	// convergence phase
	void ConvergencePhase();
	void Prediction(int frame, Frame& estPos, deque<double>& estInt);
	vector<double> Polyfit(Track tracks, string direction, int datapoints, int polydegree);		// predictor for convergence phase
	/*
	 * Function: predict the next point with Wiener Predictor using LMS algorithm
	 * Input: tracks: the track to be predicted
	 * 		  direction: to indicate the axis to be predicted
	 * 		  order: the order of the Wiener Predictor
	 * Notice: the order should be less than the length of track.
	 * Output: the next point
	 */
	double LMSWienerPredictor(Track tracks, string direction, int order);
	void Shake(Frame& estimate, deque<double>& intensity);
	deque<int> Rem(Frame& pos3D, deque<double>& int3D, double mindist_3D);
	Frame IPRonResidual(Calibration& calib, Tiff2DFinder& t, deque<int**>& pixels_orig, deque<int**>& pixels_reproj, deque<int**>& pixels_res, Frame& estimates);
	void MakeShortLinkResidual(int nextFrame, Frame& candidates, deque<Track>::iterator& tr, int iterations, bool* erase, bool* candidate_used);
	void MakeShortLinkResidual(int nextFrame, Frame& candidates, deque<Track>::iterator& tr, int iterations);

//	bool CheckVelocity(deque<Track>::iterator& tr);
	bool CheckVelocity(Track& tr);
	bool CheckAcceleration(deque<Track>::iterator& tr);
	void GetAccThred();
//	bool CheckLinearFit(deque<Track>::iterator& tr);
	bool CheckLinearFit(Track& tr);

	void MatTracksSave(string addres, string s, bool is_back_STB);
	void MatfileSave(deque<Track> tracks, string address, string name, int size);
	void SaveTrackToTXT(deque<Track> tracks, string address);
	void LoadAllTracks(string address, string frame_number, bool is_back_STB);
	void LoadTrackFromTXT(string path, TrackType trackType);

	//###################### TEMPORARY FUNTIONS FOR TESTING ###############################

	// to load 3D positions without IPR
	Frame Load_3Dpoints(string path);

	// tiff image address
	string tiffaddress;

//private:

	static const int UNLINKED = -1;
	static const int MINTRACK = 10;
	static const int EMPTY = INT_MAX;

	IPR _ipr;
	OTF OTFcalib;

	// first and last frames
	int first, last;
	deque<string> imgNameFiles;
	deque< deque<string> > imgSequence;

	// ################## initialization phase #################
	string iprfile;
	string pfieldfile;
	vector<Frame> iprMatched;
	bool iprFlag;	//use ipr or .mat files for 3D positions?
	double searchRadiusSTB; //mm 

	// ################# convergence phase ###################
	
	deque<int**> pixels_orig;									// original images
	deque<int**> pixels_reproj;									// repojected images
	deque<int**> pixels_res;									// residual image
	deque<Frame> iframes;										// 2D positons
	double avgIPDist;											// avg. interparticle distance
	double largestShift;										// Largest expected particle shift
	double fpt;													// a multiplying factor on particle intensity in order to ensure the residual has no traces of tracked particles
	double intensityLower;
	int it_innerloop;											// no. of innerloop iterations
	double maxAbsShiftChange;									// maximum allowed change in particle shift accross 2 frames
	double maxRelShiftChange;									// maximum allowed relative achange in particle shift accross 2 frames
	double acc_diff_thred = 0;
	bool enable_acc_check = false;

	// camera parameters
	std::deque<int> camID;
	std::deque<Camera> cams;
	std::deque<Camera> camsReduced;
	int ncams, Npixh, Npixw;
	bool reducedCams;

	// storing the tracks
	deque<Track> activeLongTracks;								// tracks with more than 3 particles
	deque<Track> activeShortTracks;								// tracks with 3 or less particles
	deque<Track> inactiveTracks;								// tracks that are inactive as they could not find the correct link
	deque<Track> exitTracks;									// tracks that left the measurement domain (out of at least 2 cameras)
	deque<Track> inactiveLongTracks;
	deque<Track> bufferTracks;							//  new long tracks added from Back or Forward STB (multi-pass)

	void Load_Tracks(string path, TrackType trackType);
	// dummy variables to identify the no. of tracks added and subtracted 
	int a_as = 0, a_al = 0, a_is = 0, s_as1 = 0, s_as2 = 0, s_as3 = 0, s_as4 = 0, s_al = 0, a_il = 0;
	//TESTING
	Frame tempPredictions;
	double particle_search_radius = 0;
	double designated_search_radius = 5;
	double linear_fit_error[100];
};

inline STB::STB(STB& s) {
	activeLongTracks = s.activeLongTracks;
	activeShortTracks = s.activeShortTracks;
	inactiveTracks = s.inactiveTracks;
	inactiveLongTracks = s.inactiveLongTracks;
	exitTracks = s.exitTracks;
}
#endif

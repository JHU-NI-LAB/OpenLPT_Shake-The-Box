/*
 *  Calibration.h
 *
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 */

#ifndef CALIBRATION_H
#define	CALIBRATION_H

#include <string>
#include <deque>
#include <stdexcept>
#include <utility>

#include <Camera.h>
#include <Frame.h>
#include <Position.h>
#include "Common.h"

using namespace std;
class Calibration {
public:
	// build a Calibration object from a file
	Calibration(std::string& fname, int ignoredCam, double mindist2D, double mindist3D, int nCams) throw(runtime_error);
	Calibration(deque<Camera>& cams, double mindist2D, double mindist3D, int nCams) throw(runtime_error);
	~Calibration() {
		//cout << cams[3].Get_Npixw() << endl;
		for (int n = 0; n < rcams; n++)
		{
			for (int i = 0; i < cams[n].Get_Npixh(); i++)
			{
				delete[] is_Particle[n][i];
				delete[] is_ParticlePos[n][i];
			}
			delete[] is_Particle[n];
			delete[] is_ParticlePos[n];
		}
		delete[] is_Particle;
		delete[] is_ParticlePos;
		if (!(debug_mode == SKIP_IPR_TRIANGULATION || debug_mode == SKIP_IPR_SHAKING || debug_mode == SKIP_PREVIOUS_TRACKS ||
				debug_mode == SKIP_PREVIOUS_BACK_STB)) { // when skipping tragulation, then these two variables won't be asigned
														// And thus no need to free its memory.
			delete[] camID;
			delete[] rID;
		}
		//cout << "deleted is_Particle and is_ParticlePos" << endl;
	};

	void writeGDFHeader(std::string filename);
	void fixHeader(int nr, int cols);
	// do the stereomatching
	Frame Stereomatch(const std::deque<Frame>& iframes, int framenumber, int ignoreCam) throw(std::runtime_error);

	// create an is_Particle matrix of boolean where true represents the presence of a particle in that pixel
	void Peak_identifier(int camID, int rID, Frame& f, Frame& f2);

	struct ParallelLines;
	struct CamParam;

	pair<double, double> LineParam(int camid1, int camid2, Position P);

	// functions to find the particles within mindist_2D around the projected line
	void ParticleFinder1to1(CamParam C, ParallelLines P, Position center, Position perpdir, int k, vector<Frame::const_iterator>& tmplist);
	void PixelSearchx(CamParam C, ParallelLines P, double ypix, int& Min, int& Max);
	void PixelSearchy(CamParam C, ParallelLines P, double ypix, int& Min, int& Max);
	// roots for intersection of a projected line(mm) with the pixel image 
	double Ypix1(CamParam P, double m, double c, double xpix);
	double Ypix2(CamParam P, double m, double c, double xpix);
	double Xpix1(CamParam P, double m, double c, double ypix);
	double Xpix2(CamParam P, double m, double c, double ypix);

	// function to find particles around the intersection point of 2 projected lines
	void ParticleFinder2to1(int camID, int rID, int camid1, int camid2, Position P1world, Position P2world, double mindist_1D, vector<Frame::const_iterator>& tmplist);
	// functions to check if particle candidates are actual potential candidates
	bool ParticleCheck2to1(int camid0, int camid1, int camid2, Position P1world, Position P2world, Frame::const_iterator pA, double mindist_1D);
	bool ParticleCheck1to1(int camid0, int camid1, Position P0world, Frame::const_iterator P1, double mindist_2D);

	// Loading the cleanslist from mat files
	void Load_cleanlist(string path, int frame, deque<Frame>& corrFrames_pixel, deque<Frame>& corrFrames, vector<Frame::const_iterator>*& cleanlist);

	//get cam parameters
	std::deque<Camera> Get_cam();

	//get mindist_3D and mindist_2D
	double Get_min3D();
	double Get_min2D();

	//set min3D and min 2D
	void Set_min3D(double dist);
	void Set_min2D(double dist);
	
	deque<Position> good2Dpos;
private:
	//is_Particle represents the presence of particles in a pixel and is_ParticlePos stores the particle positions in their respective pixels_orig
	bool*** is_Particle;
	vector<Frame::const_iterator>*** is_ParticlePos;

	std::string outname;
	std::ofstream outfile;
	int ncams, rcams;
	int *camID, *rID;
	std::deque<Camera> cams;
	// threshold distance between a line of sight and real particles *on each image plane* (in mm)
	double mindist_2D;
	// threshold distance between nearby lines of sight (in 3D, in mm) to match a particle
	double mindist_3D;

	// create a 3D world position from multiple positions on image planes
	std::pair<double,Position> WorldPosition(std::deque<Position> ipos, int ignoreCam) throw(std::runtime_error);

	int GroupAndPickMin( deque<double>& raydists, deque< deque<int> >& frame_index, int* buffer,
		 int num_particle, int num_match, int camera_num);

};

inline std::deque<Camera> Calibration::Get_cam() {
	return cams;
}

inline double Calibration::Get_min3D() {
	return mindist_3D;
}

inline double Calibration::Get_min2D() {
	return mindist_2D;
}

inline void Calibration::Set_min3D(double dist) {
	mindist_3D = dist;
}

inline void Calibration::Set_min2D(double dist) {
	mindist_2D = dist;
}

#endif // CALIBRATION_H

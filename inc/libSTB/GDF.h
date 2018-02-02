/*
 *  GDF.h
 *
 */

#ifndef GDF_H
#define GDF_H

#include <string>
#include <deque>
#include <stdexcept>

#include <Frame.h>

class GDF{

public:
	// constructor: process the given pixel array
	GDF(std::string filename, int start) throw(std::out_of_range);
	// destructor
	~GDF();
	
	// read GDF files
	int readGDF2D(int frame);

	// write GDF files
	void writeGDF3D(float pX, float pY, float pZ, int framenumber);
	
	// fix header information
	void fixHeader(int nr, int cols);
	
	// make a Frame object with the particle positions.
	Frame CreateFrame();
	// return the number of particles found
	int NumParticles() const;;
	
private:

	std::string outname;
	std::ifstream infile;
	std::ofstream outfile;
	double filePos[3];

	int magic,tmpi, cols, rows;
	float xi, yi, fi;
	
	int prevFrameNum;
	int currentFrameNum;
	int nextFrameNum;
	int missedFrame;
	
	int waiting_to_be_written;

	// store vectors of the x and y coordinates of the particles
	std::deque<double> x;
	std::deque<double> y;

};

inline int GDF::NumParticles() const
{
  return x.size();
}

inline GDF::~GDF()
{
  infile.close();
  outfile.close();
}

#endif // GDF_H

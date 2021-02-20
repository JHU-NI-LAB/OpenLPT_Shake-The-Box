//
//  PhotronAVI.h
//  
//  This class is a replacement for the MCine class, and should work as a 
//  seamless replacement.
//
//  Created by Nicholas Ouellette on 9/27/11.
//  Copyright 2011 Yale University. All rights reserved.
//

#ifndef WESLEYANCPV_H
#define WESLEYANCPV_H

#include <string>
#include <fstream>
#include <stdexcept>

#define BUFFERSIZE 4

enum position {
	first, current, tmp
};

class WesleyanCPV {
public:
	// constructor: give a filename
	WesleyanCPV(const std::string& name, int start = -1, int end = -1) throw(std::runtime_error, std::out_of_range);
	// destructor
	~WesleyanCPV();

	// get the parameters
	int Threshold() const;
	int Rows() const;
	int Cols() const;
	int Frames() const;
	int Colors() const;

	// get the next frame.
	int DecodeNextFrame(int** pixels, int frame) throw(std::runtime_error, std::out_of_range);

private:
	// the filename
	std::string filename;
	// the input stream itself
	std::ifstream file; 
	int threshold;
	// number of rows and columns
	short int cols;
	short int rows;
	// number of frames in the movie
	int nframes;

	int prevFrameNum;
	int currentFrameNum;
	int nextFrameNum;
	int missedFrame;
	int numPixels;
	
	int waiting_to_be_written;
	
	double filePos[3];

	char Buffer[BUFFERSIZE];
	unsigned char* buffer;

	// assume 8-bit images
	static const int DEPTH = 1;

	void Open() throw(std::runtime_error);

};

// inline functions
inline int WesleyanCPV::Threshold() const
{
  return threshold;
}

inline int WesleyanCPV::Rows() const
{
  return rows;
}

inline int WesleyanCPV::Cols() const
{
  return cols;
}

inline int WesleyanCPV::Frames() const
{
  return nframes;
}

inline int WesleyanCPV::Colors() const
{
  return ((1 << (8 * DEPTH)) - 1);
}

#endif // WESLEYANCPV_H

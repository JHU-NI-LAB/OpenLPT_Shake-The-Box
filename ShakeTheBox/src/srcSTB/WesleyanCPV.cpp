//
//  WesleyanCPV.cpp
//  
//
//  Created by Nicholas Ouellette on 9/27/11.
//  Copyright 2011 Yale University. All rights reserved.
//  
//	Modified by Stefan Kramel
//	

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <string.h>

#include <WesleyanCPV.h>

using namespace std;

WesleyanCPV::WesleyanCPV(const string& name, int start, int end) throw(runtime_error, out_of_range)
: filename(name)
{
	try {
		Open();
	} 
	catch (runtime_error& e) {
		cerr << e.what() << endl;
		throw runtime_error("runtime_error caught from WesleyanCPV::Open()");
	}
	
	// if desired, seek to the starting frame
	filePos[first] = file.tellg();
	file.read(reinterpret_cast<char*>(&currentFrameNum), 4);
	file.read(reinterpret_cast<char*>(&Buffer), 4);
	numPixels = ((unsigned char)Buffer[3] << 14) + ((unsigned char)Buffer[2] << 6) + ((unsigned char)Buffer[1] >> 2);
	file.seekg((numPixels*4), ios::cur);
	
	while (currentFrameNum != start) {
		filePos[first] = file.tellg();
		file.read(reinterpret_cast<char*>(&currentFrameNum), 4);
		if (currentFrameNum == 0) // end of frame
			break;
		file.read(reinterpret_cast<char*>(&Buffer), 4);
		numPixels = ((unsigned char)Buffer[3] << 14) + ((unsigned char)Buffer[2] << 6) + ((unsigned char)Buffer[1] >> 2);
		file.seekg((numPixels*4), ios::cur);
	}
	file.seekg(filePos[first], ios::beg);
	
	waiting_to_be_written = 0;
	
	nframes = end - start;
}

WesleyanCPV::~WesleyanCPV()
{
  file.close();
  delete []buffer;
}

int WesleyanCPV::DecodeNextFrame(int** pixels, int frame) throw(runtime_error, out_of_range)
{	
	
	while(!waiting_to_be_written) {
		filePos[current] = file.tellg();
		file.read(reinterpret_cast<char*>(&currentFrameNum), 4);
		cout << "currentFrameNum: " << currentFrameNum << endl;
		file.read(reinterpret_cast<char*>(&Buffer), 4);
		numPixels = ((unsigned char)Buffer[3] << 14) + ((unsigned char)Buffer[2] << 6) + ((unsigned char)Buffer[1] >> 2);
		file.seekg((numPixels*4), ios::cur);
		filePos[tmp] = file.tellg();
		file.read(reinterpret_cast<char*>(&nextFrameNum), 4);
	
		cout << "nextFrameNum: " << nextFrameNum << endl;
	
		file.seekg(filePos[current], ios::beg);	

		if (prevFrameNum == currentFrameNum-1 || nextFrameNum == currentFrameNum+1) {
			file.read(reinterpret_cast<char*>(&currentFrameNum), 4);
			file.read(reinterpret_cast<char*>(&Buffer), 4);
			numPixels=((unsigned char)Buffer[3]<<14)+((unsigned char)Buffer[2]<<6)+((unsigned char)Buffer[1]>>2);
			for (int i = 0; i < numPixels; i++) {	
				file.read(reinterpret_cast<char*>(&Buffer), 4);
				buffer[cols*(((unsigned char)Buffer[2]>>3)+(((unsigned char)Buffer[3])<<5))+((unsigned char)Buffer[1]+(((unsigned char)Buffer[2]&07)<<8))]=(unsigned char)Buffer[0];
			}
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					pixels[i][j] = static_cast<int>(buffer[cols * i + j]);
					memset(&buffer[cols * i + j], 0, sizeof(unsigned char));
				}
			}
			waiting_to_be_written = 1;
			prevFrameNum = currentFrameNum;
			break;
		}
		else {
			file.seekg(filePos[tmp], ios::beg);
		}
	}
	if(currentFrameNum == frame) {
		missedFrame = 0;
		waiting_to_be_written = 0;
	}
	else {
		missedFrame = 1;
	}
		
	return (missedFrame);
}

// Open .cpv file, based on CPVPlayer decoder.cpp
void WesleyanCPV::Open() throw(runtime_error)
{
	try {
		file.open(filename.c_str(), ios::in | ios::binary);
	} 
	catch (...) {
		throw runtime_error("Cannot open .CPV file!");
	}
	if (file.good()) {
		file.seekg(0, ios::beg); //seek to first bit
		file.seekg(4, ios::cur); //skip reading version & numlines (1 byte each)
		file.read(reinterpret_cast<char*>(&cols), 2);
		file.read(reinterpret_cast<char*>(&rows), 2);
		file.seekg(12, ios::cur); //skip reading exptime, fps & gain (4 byte each) we are at position (20, ios::beg)
	}
	else {
		cerr << "Failed to open file" << endl;
		exit (1);
	}
	// finally, allocate the pixel buffer
	try {
		buffer = new unsigned char[rows * cols];
		
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				memset(&buffer[cols * i + j], 0, sizeof(unsigned char));
				if (static_cast<int>(buffer[cols * i + j]) != 0){
					cerr << "Failure to empty buffer in WesleyanCPV::Open" << endl;
					exit(1);
				}
			}
		}
		
	} 
	catch (bad_alloc&) {
		throw runtime_error("Failure to allocate buffer in WesleyanCPV::Open");
	}	
}

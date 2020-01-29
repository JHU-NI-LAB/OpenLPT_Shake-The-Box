/*
 *  GDF.cpp
 */

#include <iostream>
#include <fstream>
#include <queue>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <GDF.h>
#include <Position.h>
#include <WesleyanCPV.h>

using namespace std;

GDF::GDF(std::string filename, int start) throw(out_of_range)
: outname(filename)
{
	// Read Header:
	infile.open(outname.c_str(), ios::in | ios::binary);
	if (infile.is_open()){
		infile.read(reinterpret_cast<char*>(&magic), 4);
		// number of dimensions
		infile.read(reinterpret_cast<char*>(&tmpi), 4);
		// number of columns
		infile.read(reinterpret_cast<char*>(&cols), 4);
		// number of rows
		infile.read(reinterpret_cast<char*>(&rows), 4);
		// 4 means floating point numbers (see Matlab read_gdf function for more info)
		infile.read(reinterpret_cast<char*>(&tmpi), 4);
		// number of total points
		infile.read(reinterpret_cast<char*>(&tmpi), 4);
	}
	filePos[first] = infile.tellg();
	infile.seekg(20, ios::cur);
	infile.read(reinterpret_cast<char*>(&fi), 4);
	currentFrameNum = fi;
	
	while(currentFrameNum != start) {
		filePos[first] = infile.tellg();
		infile.seekg(20, ios::cur);
		infile.read(reinterpret_cast<char*>(&fi), 4);
		currentFrameNum = fi;
	}
		
	//cout << "currentFrameNum: " << currentFrameNum << endl;	
	
	infile.seekg(filePos[first], ios::beg);

	waiting_to_be_written = 0;

}

int GDF::readGDF2D(int frame) {
	//deque<double> newx;
	//deque<double> newy;
	
    filePos[current] = infile.tellg();
    infile.seekg(20, ios::cur);
    infile.read(reinterpret_cast<char*>(&fi), 4);
    currentFrameNum = fi;
    
    if(currentFrameNum == frame) {
        missedFrame = 0;
        waiting_to_be_written = 0;
    }
    else {
        missedFrame = 1;
        infile.seekg(filePos[tmp], ios::beg);
        return (missedFrame);
    }
    
    x.clear();
    y.clear();

    
    infile.seekg(filePos[current], ios::beg);
    
    int Nparticle=0;
    
    while (fi==currentFrameNum) {
        
        //cout << "\tParticlecount: " << i << endl;
        infile.read(reinterpret_cast<char*>(&xi), 4);
        x.push_back(xi);
        infile.read(reinterpret_cast<char*>(&yi), 4);
        y.push_back(yi);

        infile.seekg(12, ios::cur); //skip 12 bytes
        infile.read(reinterpret_cast<char*>(&fi), 4);
        if (infile.eof()) {
            std::cout << "\tEnd of file reached." << endl;
            break;
        }
        Nparticle++;
        //x = newx;
        //y = newy;
    }
    
    cout << "first particle:" << x.at(0) << ",\t" << y.at(0) << "\n";

    
    if (fi>currentFrameNum) {
        infile.seekg(-24, ios::cur); //skip 16 bytes
        x.pop_back();
        y.pop_back();


    }
  
    return (missedFrame);

    
    
	/*
     while (!waiting_to_be_written) {
		int particlecount = 0;
		filePos[current] = infile.tellg();
		infile.seekg(20, ios::cur);
		infile.read(reinterpret_cast<char*>(&fi), 4);
		currentFrameNum = fi;
        
        
		//cout << "prevFrameNum: " << prevFrameNum << endl;
		//cout << "currentFrameNum: " << currentFrameNum << endl;
		infile.seekg(filePos[current], ios::beg);
        
        
		while(currentFrameNum == fi){
			particlecount++;
			filePos[tmp] = infile.tellg();
			infile.seekg(20, ios::cur);
			infile.read(reinterpret_cast<char*>(&fi), 4);
			nextFrameNum = fi;
			if (infile.eof()) {
				std::cout << "\tEnd of file reached." << endl;
				break;
			}
		}
		particlecount -= 1;
		//cout << "nextFrameNum: " << nextFrameNum << endl;
		//cout << "\tParticlecount: " << particlecount << endl;
		infile.seekg(filePos[current], ios::beg);
		
		if (prevFrameNum == currentFrameNum-1 || nextFrameNum == currentFrameNum+1) {
			for(int i = 0; i < particlecount; i++){
                //cout << "\tParticlecount: " << i << endl;
				infile.read(reinterpret_cast<char*>(&xi), 4);
				newx.push_back(xi);
				infile.read(reinterpret_cast<char*>(&yi), 4);
				newy.push_back(yi);
				infile.seekg(16, ios::cur); //skip 16 bytes
				x = newx;
				y = newy;
			}
		waiting_to_be_written = 1;
		prevFrameNum = currentFrameNum;
		break;
		}
		else {
			infile.seekg(filePos[tmp], ios::beg);
		}
     
	}*/
    

}

void GDF::fixHeader(int nr, int cols){
	
	// now fix up the header with the proper sizes
	outfile.seekp(8, ios::beg);
	outfile.write(reinterpret_cast<const char*>(&cols), 4);
	outfile.write(reinterpret_cast<const char*>(&nr), 4);
	outfile.seekp(4, ios::cur);
	int tmpi = cols * nr;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);

}

Frame GDF::CreateFrame() {
  deque<Position> pos;
  deque<double>::const_iterator xi_end = x.end();
  deque<double>::const_iterator yi = y.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
    pos.push_back(Position(*xi, *yi, 0));
  }
	return Frame(pos);
}

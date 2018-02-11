/*
* Based on STB (Schanz et. al) tracking code,
*
* Written 2/28/17 by Ashwanth
*
* Latest version: 9/28/2017
*
*/

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <deque>
#include <vector>

#include <GDF.h>
#include <Frame.h>
#include <STB.h>
#include <gnu/libc-version.h>


using namespace std;

// configuration parameters
struct ConfigFile {
	int ncams;
	deque<int> camIDs;
	deque<string> imgNameFiles;
	string iprfile;
	string pfieldfile;
	//double fps;
	double threshold;
	//double cluster_rad;
	//int npredict;
	//double max_disp;
	//int memory;
	int first;
	int last;
	string stereomatched;
	string outname;
	bool iprFlag;
	double initialPhaseRadius;
	double avgSpace;
	double largestShift;
	double maxAbsShiftChange;
	double maxRelShiftChange;
	double fpt;
	double lowerInt;
};

// globals
struct ConfigFile config;

void ImportConfiguration(struct ConfigFile* config, char* name);

 int main(int argc, char** argv) {
//	 printf("GNU libc version: %s\n", gnu_get_libc_version());
		if (argc < 2) {

			cerr << "Usage: " << argv[0] << " <configuration file>" << endl;
		exit(1);
	}

	ImportConfiguration(&config, argv[1]);

	// read the camera calibration information
	//Calibration calib(config.iprfile);
	//int first = config.first;
	//int last = config.last;
	//int threshold = config.threshold;
	/*double cluster_rad = config.cluster_rad;*/

	// do the stereomatching
	//cout << "Stereomatching..." << endl;
	//calib.writeGDFHeader(config.stereomatched);


	// tracking using STB
	STB s(config.first, config.last, config.pfieldfile, config.iprfile, config.ncams, config.camIDs, config.imgNameFiles,
		config.initialPhaseRadius, config.avgSpace, config.largestShift, config.maxAbsShiftChange,
		config.maxRelShiftChange, config.fpt, config.lowerInt, config.iprFlag);
	// saving the tracks
	s.MatTracksSave(s.tiffaddress,"", config.last);
	// applying a pass of BackSTB
/*	BackSTB bs(first, last, config.initialPhaseRadius);
	bs.UpdateTracks(s);
	// saving the tracks
	s.MatTracksSave("back");
	// applying a pass of ForwardSTB
	ForwardSTB fs(first, last, config.initialPhaseRadius);
	fs.UpdateTracks(s);
	// saving the tracks
	s.MatTracksSave("forward");
*/
	//std::cout << "\tTotal number of stereomatched particles: " << nr << endl;
	// first argument: total number of particles; second argument: number of columns in .gdf file;
	// framenumber, x, y, z, intersect, xy, xy, xy, xy;
	//calib.fixHeader(nr,13);

	// Output 3D positions in .txt file
	/*
	// 3D Position output for each camera, saved in pos3Dcam[camID].txt
	std::ostringstream fn1;
	fn1 << "pos3Dold.txt";
	std::ofstream outfile1;
	outfile1.open (fn1.str().c_str(),std::ios_base::binary);
	for (int i = 0; i < (config.last - config.first); ++i){
		outfile1 << "Frame: " << i << endl;
		outfile1 << matched[i] << endl;
	}
	outfile1.close();
	*/



	// finally, do the tracking
/*	cout << "Tracking..." << endl;
	Tracker::TrackMode mode = Tracker::FRAME3;

	if (config.npredict == 0) {
		mode = Tracker::FRAME2;
	}
	else if (config.npredict == 1) {
		mode = Tracker::FRAME3;
	}
	else if (config.npredict == 2) {
		mode = Tracker::FRAME4;
	}

	Tracker t(mode, config.max_disp, config.memory, config.fps,
	config.outname);
	t.MakeTracks(matched);*/
	// Done!
	cout << "Done." << endl;
	return 0;
}

void ImportConfiguration(struct ConfigFile* config, char* name) {
		cout << "Reading configuration file..." << endl;
		ifstream file(name, ios::in);
		string line;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->ncams = atoi(line.c_str());

		for (int i = 0; i < config->ncams; i++) {
			getline(file, line);
			line.erase(line.find_first_of(' '));
			config->camIDs.push_back(atoi(line.c_str()));
		}

		for (int i = 0; i < config->ncams; i++) {
			getline(file, line);
			line.erase(line.find_first_of(' '));
			config->imgNameFiles.push_back(line.c_str());
		}

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->iprfile = line;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->pfieldfile = line;

		//getline(file, line);
		//line.erase(line.find_first_of(' '));
		//config->fps = atof(line.c_str());

		//getline(file, line);
		//line.erase(line.find_first_of(' '));
		//config->threshold = atof(line.c_str());

		//getline(file, line);
		//line.erase(line.find_first_of(' '));
		//config->cluster_rad = atof(line.c_str());

		//getline(file, line);
		//line.erase(line.find_first_of(' '));
		//config->npredict = atoi(line.c_str());

		//getline(file, line);
		//line.erase(line.find_first_of(' '));
		//config->max_disp = atof(line.c_str());

		//getline(file, line);
		//line.erase(line.find_first_of(' '));
		//config->memory = atoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->first = atoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->last = atoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->stereomatched = line;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->outname = line;

		getline(file, line);
		line.erase(line.find_first_of(' '));

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->iprFlag = stoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->initialPhaseRadius = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->avgSpace = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->largestShift = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->maxAbsShiftChange = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->maxRelShiftChange = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->fpt = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->lowerInt = stof(line.c_str());
}

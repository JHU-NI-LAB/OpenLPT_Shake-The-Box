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
#include <vector>
#include <ratio>
#include <chrono>
#include <sys/stat.h>
#include <GDF.h>
#include <Frame.h>
#include <STB.h>
#include "BackSTB.h"
#ifdef WINDOWS
	#include <dirent.h>
#else
	#include <gnu/libc-version.h>
#endif

#include "Common.h"
#include "BoundaryCheck.h"

using namespace std;

#ifdef WINDOWS
	char* version = "W2.1.041619"; //Version of this project
#else
	char* version = "L2.1.041619"; //Version of this project
#endif

// globals
ConfigFile config;
// debug
DebugMode debug_mode;
int debug_frame_number;
ERROR error = ERROR(0);
bool to_save_data;
// Boundary check
BoundaryCheck boundary_check;


void ImportConfiguration(struct ConfigFile* config, char* name);

void GetDebugMode() {
	cout<<"Select debug mode:\n"
			<<"0. No debug\n"
			<<"1. Load the 2D position \n"
			<<"2. Load the 3D position before shaking\n"
			<<"3. Load the 3D position after shaking\n"
			<<"4. Load the 3D positions for predictive field \n"
			<<"5. Load the predictive field \n"
			<<"6. Load the tracks from the initial phase \n"
			<<"7. Load the tracks from specific frames \n"
			<<"8. Load the tracks from specific frames for back STB \n"
			<<"Enter the NO. of the option:";
	int NO = 0;
	cin>>NO;
	if (NO < 9) debug_mode = DebugMode(NO); else debug_mode = DebugMode(0);
//	cout<<debug_mode<<endl;
	if (!(NO == 0)) {
		cout<<"Enter the frame number to be debugged:";
		cin>>debug_frame_number;
//		cout<<debug_frame_number;
	}
	cout<<"To save data(1):";
	cin>>to_save_data;
}



 int main(int argc, char** argv) {
	 printf("Code version: %s\n", version);
		if (argc < 2) {

			cerr << "Usage: " << argv[0] << " <configuration file>" << endl;
		exit(1);
	}

	ImportConfiguration(&config, argv[1]);

	boundary_check.SetLimit(config.x_upper_limt, config.x_lower_limit, config.y_upper_limt, config.y_lower_limit,
			config.z_upper_limt, config.z_lower_limit);

	GetDebugMode();

	//Create folder to save tracks
	struct stat info;

	string folder_path = config.iprfile;
	folder_path.erase(folder_path.size() - 13, 13); //erase iprconfig.txt
	folder_path = folder_path + "Tracks";

#ifdef WINDOWS
    if( stat( folder_path.c_str(), &info ) != 0 ) {
        printf( "cannot access %s\n", folder_path.c_str() );
        mkdir(folder_path.c_str());
        mkdir((folder_path + "/InitialTracks").c_str());
        mkdir((folder_path + "/ConvergedTracks").c_str());
        mkdir((folder_path + "/BackSTBTracks").c_str());
        printf( "%s have been created!\n", folder_path.c_str() );
    }
    else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
        printf( "%s exists\n", folder_path.c_str() );
    else {
            mkdir(folder_path.c_str());
            mkdir((folder_path + "/InitialTracks").c_str());
            mkdir((folder_path + "/ConvergedTracks").c_str());
            mkdir((folder_path + "/BackSTBTracks").c_str());
            printf( "%s have been created!\n", folder_path.c_str() );
    }

#else
	if( stat( folder_path.c_str(), &info ) != 0 ) {
	    printf( "cannot access %s\n", folder_path.c_str() );
	    mkdir(folder_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	    mkdir((folder_path + "/InitialTracks").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	    mkdir((folder_path + "/ConvergedTracks").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	    mkdir((folder_path + "/BackSTBTracks").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	    mkdir((folder_path + "/IPRcandidates").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	    printf( "%s have been created!\n", folder_path.c_str() );
	}
	else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
	    printf( "%s exists\n", folder_path.c_str() );
	else {
		mkdir(folder_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		mkdir((folder_path + "/InitialTracks").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		mkdir((folder_path + "/ConvergedTracks").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		mkdir((folder_path + "/BackSTBTracks").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		printf( "%s have been created!\n", folder_path.c_str() );
	}
#endif

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
//	auto start = std::chrono::system_clock::now();
	STB s(config.first, config.last, config.pfieldfile, config.iprfile, config.ncams, config.camIDs, config.imgNameFiles,
		config.initialPhaseRadius, config.avgSpace, config.largestShift, config.maxAbsShiftChange,
		config.maxRelShiftChange, config.fpt, config.lowerInt, config.iprFlag);
//	auto end = std::chrono::system_clock::now();
//	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//	cout << "Total time: "<<elapsed.count() << '\n';

	// saving the tracks
	//s.MatTracksSave(s.tiffaddress,"", config.last);
	// applying a pass of BackSTB
	if (config.backSTB) {
		BackSTB bs(config.first, config.last, config.dist_two_tracks);
		bs.UpdateTracks(s);
	}
	// saving the tracks
//	s.MatTracksSave("back");
	// applying a pass of ForwardSTB
//	ForwardSTB fs(first, last, config.initialPhaseRadius);
//	fs.UpdateTracks(s);
//	// saving the tracks
//	s.MatTracksSave("forward");

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
	cout << "All of images have been processed!" << endl;
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
		config->x_lower_limit = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->x_upper_limt = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->y_lower_limit = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->y_upper_limt = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->z_lower_limit = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->z_upper_limt = atof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->iprFlag = stoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		//unit conversion from voxel to mm
		config->factor = (config->x_upper_limt - config->x_lower_limit) / 1000;
		config->initialPhaseRadius = stof(line.c_str()) * config->factor;

		getline(file, line);
		line.erase(line.find_first_of(' '));

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->shaking_shift = stof(line.c_str()) * config->factor;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->avgSpace = stof(line.c_str()) * config->factor;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->largestShift = stof(line.c_str()) * config->factor;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->maxAbsShiftChange = stof(line.c_str()) * config->factor;

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->maxRelShiftChange = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->fpt = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->lowerInt = stof(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->backSTB = stoi(line.c_str());

		getline(file, line);
		line.erase(line.find_first_of(' '));
		config->dist_two_tracks = stof(line.c_str()) * config->factor;
}



/*
 * Common.h
 *	A common file to declare variables which are globally used
 *  Created on: Feb 26, 2018
 *      Author: shiyongtan
 */

#ifndef INC_LIBSTB_COMMON_H_
#define INC_LIBSTB_COMMON_H_
#include <string>
#include <deque>
#include "BoundaryCheck.h"
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
	double x_upper_limt, x_lower_limit;
	double y_upper_limt, y_lower_limit;
	double z_upper_limt, z_lower_limit;
	double factor; // using to converse voxel into mm.
	bool iprFlag;
	double initialPhaseRadius;
	double shaking_shift;
	double avgSpace;
	double largestShift;
	double maxAbsShiftChange;
	double maxRelShiftChange;
	double fpt;
	double lowerInt;
	bool backSTB;
	double dist_two_tracks;
};

// Data for debug
enum DebugMode {
		NO_SKIP = 0,
		SKIP_IPR_2D_POSITION,
		SKIP_IPR_TRIANGULATION,
		SKIP_IPR_SHAKING,
		SKIP_IPR,
		SKIP_PREDITIVE_FIELD,
		SKIP_INITIAL_PHASE,
		SKIP_PREVIOUS_TRACKS,
		SKIP_PREVIOUS_BACK_STB
};
extern DebugMode debug_mode;
extern int debug_frame_number;
extern bool to_save_data;
extern BoundaryCheck boundary_check;
extern ConfigFile config;

//Data for error
enum ERROR {
	NONE = 0,
	NO_FILE, //Can't find required file
};

extern ERROR error;


#endif /* INC_LIBSTB_COMMON_H_ */

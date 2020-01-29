#ifndef PredictiveField_H
#define	PredictiveField_H

#include <string>
#include <deque>
#include <stdexcept>
#include <utility>

//#include <Position.h>
#include <Frame.h>
// TODO: Temporary modification by Shiyong Tan, 2/1/18
//#include <matio.h>
//#include <OTF.h>
#include "Common.h"




using namespace std;
class PredictiveField {
public:
	// constructor : To get predictive field( takes the previous frame 3D pos, field parameters from file, frame, matlabFlag )
	PredictiveField(){};

	void GetPredictiveField(Frame prevFramePos, Frame currFramePos, std::string& fname, int frame);

	// destructor
	~PredictiveField() {
		for (int n = 0; n < 3; n++) {
			delete[] field[n];
		}
		delete[] gridPoints;
	};

	void GridPoints();
	void Field();
	void SearchVolumeParticles(deque<Frame::const_iterator>& Frame, double rsqr, vector<double>* gridPoint, int grid, int frame);
	void DisplacementMap(double*** &dispMap, deque< deque<double> > displacements);
	deque<double> DisplacementMapPeak(double*** &dispMap);
	deque<int> IndexofLargestElement(double*** &array);
	double Gaussian1DPeak(double y1, double v1, double y2, double v2, double y3, double v3);

	// getting the necessary variables
	vector< vector<double> > GetGrid() {
		vector< vector<double> > Grids;
		for (int n = 0; n < 3; n++) {
			vector<double> temp(totalGridPoints);
			for (int i = 0; i < totalGridPoints; i++) {
				temp[i] = gridPoints[n][i];
			}
			Grids.push_back(temp);
		}
		return Grids;
	}
	vector< vector<double> > GetField() {
		vector< vector<double> > Field;
		for (int n = 0; n < 3; n++) {
			vector<double> temp(totalGridPoints);
			for (int i = 0; i < totalGridPoints; i++) {
				temp[i] = field[n][i];
			}
			Field.push_back(temp);
		}
		return Field;
	}
	Frame GetCurrPos3D() {
		return matchedCurr;
	}
	Frame GetPrevPos3D() {
		return matchedPrev;
	}

	vector<double> linspace(double first, double last, int len);
	Position ParticleInterpolation(Position pos3D);

	void MatfileSave(vector<double*> pos, string name);
	void MatfileSave(vector<double> pos[3], string name);
	void MatfileSave(vector< vector<double> > pos, string name);
	void SaveField(string file_path);
	void Load_field(string path);

private:
	bool getPredictiveField;
	int currFrame, prevFrame;
	// X,Y,Z limits of view area
	double viewAreaLimits[6];
	// grid sppacing in X,Y,Z
	double gridSize[3];
	// # of grids in X,Y,Z
	int numGrids[3];
	// radius of interrogation sphere
	double radius;

	// 3D particles at prev and curr frame
	Frame matchedCurr;
	Frame matchedPrev;
	
	// grid points (X,Y,Z)
	vector<double>* gridPoints;
	int totalGridPoints;

	// parameters for displacement map
	int Size;
	double m, c;

	vector<double*> field;
	vector< vector<double> > particleDisplacements;


	// variables for data from .mat file
};



#endif

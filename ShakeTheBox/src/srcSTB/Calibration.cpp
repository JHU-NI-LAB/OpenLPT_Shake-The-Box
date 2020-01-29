/*
 *  Calibration.cpp
 *
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 */

#include <sstream>
#include <iostream>
#include <algorithm>
#include <list>
#include <vector>
#include <cmath>
#include <ctime>
#include <omp.h>
/*
 * Modified by Shiyong Tan, 2/5/18
 * Matio library has been descared. Use DataIO instead
 * Start:
 */
//#include <matio.h>
#include "NumDataIO.h"
// End

#include <Calibration.h>
#include <GDF.h>

#include "Common.h"

#define PI 3.14159
using namespace std;

// A structure to obtain cam parameters
struct Calibration::CamParam
{
	CamParam(int Npixw, int Npixh, double kr, double wpix, double hpix, double Noffh, double Noffw)
		: Npixw(Npixw), Npixh(Npixh), kr(kr), wpix(wpix), hpix(hpix), Noffh(Noffh), Noffw(Noffw)
	{}
	// camera parameters
	int Npixw, Npixh;
	double kr, wpix, hpix, Noffh, Noffw, c_proj, m_proj, mindist_2D;
};

// A structure to define the two lines at +-mindist_2D distance around the projected line
struct Calibration::ParallelLines
{
	ParallelLines(double C, double m, double del)
		: c_proj(C), m_proj(m), mindist_2D(del)
	{}
	double c_proj, m_proj, mindist_2D;
	
		// theta1 is the angle perpendicular to the projected line
	double theta1 = atan(-1 / m_proj);
			
		// the intercepts of parallel lines y = mx+c (lines at +-mindist_2D distance from proj line)
	double c_plusdel = c_proj + mindist_2D*(sin(theta1)-m_proj*cos(theta1));
	double c_minusdel = c_proj - mindist_2D*(sin(theta1) - m_proj*cos(theta1));
};

// Constructor
Calibration::Calibration(std::string& fname, int ignoredCam, double mindist2D, double mindist3D, int nCams) throw(runtime_error) 
	: mindist_2D(mindist2D), mindist_3D(mindist3D), ncams(nCams)
{
	// remove comments from the file
	ifstream infile(fname.c_str(), ios::in);
	string line;
	stringstream parsed;
	while (getline(infile, line)) {
		size_t commentpos = line.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos) {
				line.erase(commentpos);
			}
			parsed << line << '\t';
		}
	}
	infile.close();

	// read the number of cameras
	parsed >> nCams;
	if (ignoredCam < ncams) { rcams = ncams - 1; }
	else					 { rcams = ncams; }

	// now read in the parameters for each camera
	for (int i = 0; i < ncams; ++i) {
		cams.push_back(Camera(parsed));
	}

	camID = new int[rcams]; rID = new int[rcams];
	// initializing the logic matrix (is_Particle) and its corresponding particle vectors (is_ParticlePos) 
	Calibration::is_Particle = new bool**[rcams];
	Calibration::is_ParticlePos = new vector<Frame::const_iterator>**[rcams];

	for (int n = 0; n < rcams; n++)	{
		int Npixw = cams[n].Get_Npixw();
		int Npixh = cams[n].Get_Npixh();
		Calibration::is_Particle[n] = new bool*[Npixh];
		try {
			Calibration::is_ParticlePos[n] = new vector<Frame::const_iterator>*[Npixh];
		}
		catch (bad_alloc&) {
			runtime_error("cannot assign storage space for is_ParticlePos");
		}
		for (int i = 0; i < Npixh; i++)	{
			Calibration::is_Particle[n][i] = new bool[Npixw];
			try {
				Calibration::is_ParticlePos[n][i] = new vector<Frame::const_iterator>[Npixw];
			}
			catch (bad_alloc&) {
				runtime_error("cannot assign storage space for is_ParticlePos");
			}
		}
	}
}

Calibration::Calibration(deque<Camera>& cams, double mindist2D, double mindist3D, int nCams) throw(runtime_error)
	: cams(cams), mindist_2D(mindist2D), mindist_3D(mindist3D), ncams(nCams)
{
	rcams = ncams;

	camID = new int[rcams]; rID = new int[rcams];

	// initializing the logic matrix (is_Particle) and its corresponding particle vectors (is_ParticlePos) 
	Calibration::is_Particle = new bool**[rcams];
	Calibration::is_ParticlePos = new vector<Frame::const_iterator>**[rcams];

	for (int n = 0; n < rcams; n++) {
		int Npixw = cams[n].Get_Npixw();
		int Npixh = cams[n].Get_Npixh();
		Calibration::is_Particle[n] = new bool*[Npixh];
		try {
			Calibration::is_ParticlePos[n] = new vector<Frame::const_iterator>*[Npixh];
		}
		catch (bad_alloc&) {
			runtime_error("cannot assign storage space for is_ParticlePos");
		}
		for (int i = 0; i < Npixh; i++) {
			Calibration::is_Particle[n][i] = new bool[Npixw];
			try {
				Calibration::is_ParticlePos[n][i] = new vector<Frame::const_iterator>[Npixw];
			}
			catch (bad_alloc&) {
				runtime_error("cannot assign storage space for is_ParticlePos");
			}
		}
	}
}

void Calibration::writeGDFHeader(std::string outname)
{
	// open the output file with a temporary header
	outfile.open(outname.c_str(), ios::out | ios::binary);
	int magic = 82991;
	outfile.write(reinterpret_cast<const char*>(&magic), 4);
	// number of dimensions
	int tmpi = 2;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// number of columns
	tmpi = 6;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// number of rows: we don't know this yet!
	tmpi = 0;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// a 4 means floating point numbers
	tmpi = 4;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	// number of total points: we don't know this yet!
	tmpi = 0;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);

	//std::cout << "\nHeader information written..." << endl;

}

void Calibration::fixHeader(int nr, int cols){

	// now fix up the header with the proper sizes
	outfile.seekp(8, ios::beg);
	outfile.write(reinterpret_cast<const char*>(&cols), 4);
	outfile.write(reinterpret_cast<const char*>(&nr), 4);
	outfile.seekp(4, ios::cur);
	int tmpi = cols * nr;
	outfile.write(reinterpret_cast<const char*>(&tmpi), 4);
	//std::cout << "\nHeader information updated!" << endl;
}

/*
 * ********************************************************************************************************
 * Following is the merge sort algorithm designed for ranging the matches according to the match preference
 */
void CopyArray(int* A, int iBegin, int iEnd, int* B) {
	for (int k = iBegin; k < iEnd; ++k) {
		B[k] = A[k];
	}
}

void Merge(int* A, int iBegin, int iMiddle, int iEnd, int* B, double** C) {
	int i = iBegin, j = iMiddle;
	for(int k = iBegin; k < iEnd; ++k) {
		if (i < iMiddle && j >= iEnd) {
			B[k] = A[i];
			i = i + 1;
		} else if (i < iMiddle && C[A[i]][0] >= C[A[j]][0]) {
			if (C[A[i]][0] > C[A[j]][0]) {
				B[k] = A[i];
				i = i + 1;
			} else { // when the preference numbers are the same, then compare the error
				if (C[A[i]][1] < C[A[j]][1]) {  // the smaller error should rank the first
					B[k] = A[i];
					i = i + 1;
				} else {
					B[k] = A[j];
					j = j + 1;
				}
			}
		} else {
			B[k] = A[j];
			j = j + 1;
		}
	}
}

void SplitMerge(int* B, int iBegin, int iEnd, int* A, double** C) {
	if (iEnd - iBegin < 2) return;
	int iMiddle = (iEnd + iBegin) / 2;
	SplitMerge(A, iBegin, iMiddle, B, C);
	SplitMerge(A, iMiddle, iEnd, B, C);
	Merge(B, iBegin, iMiddle, iEnd, A, C);
}

void MergeSort(int* A, double** C, int n) {
	int B[n];
	CopyArray(A, 0, n, B);
	SplitMerge(B, 0, n, A, C);
}
//********************************************* Merge Sort Ends here ***************************************************//

Frame Calibration::Stereomatch(const deque<Frame>& iframes, int framenumber, int ignoreCam = 100) throw(runtime_error) {

	vector < deque<Position> > PosTouse;
/*	if (iframes.size() != cams.size()) {
		throw runtime_error("Number of cameras and number of images do not match!");
	} */

	// correct all the distortions, and move the points into a coordinate
	// system (in mm) with the origin in the middle of the image plane.
	deque<Frame> corrframes; // physical units
	deque<Frame> corrframes_Pixel; // pixel units
	double mindist_1D = mindist_2D;
	int temps = 0;
	//std::cout << "\tCorrecting distortion..." << endl;


	int id = 0;
	for (int n = 0; n < ncams; n++) {
		if (n != ignoreCam) {
			camID[id] = n;
			id++;
		}
	}
	// for reduced cams
	for (int n = 0; n < rcams; n++) {
			rID[n] = n;
	}

	for (int i = 0; i < rcams; ++i) {
		deque<Position> corrpos;
		deque<Position> corrpos_Pixel;
		Frame::const_iterator fitend = iframes[i].end();
		//int n = 0;
		for (Frame::const_iterator fit = iframes[i].begin(); fit != fitend; ++fit) {
			/*
			 * Modified by Shiyong Tan
			 * date: 6.8.18
			 * Cams[i] is not correct for reduced camera
			 * Start:
			 */
//			corrpos.push_back(cams[i].UnDistort(*fit));
			corrpos.push_back(cams[camID[i]].UnDistort(*fit));
			// End
			corrpos_Pixel.push_back(*fit);
			//n = n + 1;
		}
		corrframes.push_back(Frame(corrpos));
		corrframes_Pixel.push_back(Frame(corrpos_Pixel));
	}

//	double duration;
//	clock_t start;
//	clock_t start1;
//	start = clock();

//	camID = new int[rcams]; rID = new int[rcams];
//	int id = 0;
//	for (int n = 0; n < ncams; n++) {
//		if (n != ignoreCam) {
//			camID[id] = n;
//			id++;
//		}
//	}
//	// for reduced cams
//	for (int n = 0; n < rcams; n++) {
//			rID[n] = n;
//	}

	// creating the logic matrix and its corresponding particle lists
	for (int i = 0; i < rcams; i++)	
		Peak_identifier(camID[i], rID[i], corrframes_Pixel[i], corrframes[i]);
	

	// for 1st camera, draw a line of sight through each particle on its image plane.
	// project these lines of sight onto the image planes of 2nd camera.
	// particles within mindist_2D of these lines are candidate matches between first 2 cams.
	// then project 2 line of sights from each particle pair of 1st n 2nd cam onto 3rd cam.
	// particles within mindist_1D of the intersection point are candidate matches from 3rd cam.
	// repeat similarly for subsequent cams
	//std::cout << "\tConstructing clean list..." << endl;
	int avgsize = 0;
	int numlists = 0;
	
	// a cleanlist that contains matched particle list from each camera (up to 4 cameras)
	int nCams = (rcams < 4) ? rcams : 4;
	vector<Frame::const_iterator>* cleanlist = new vector<Frame::const_iterator> [nCams];

	// looping over each particle in 1st camera to find a match in subsequent cameras
	Frame::const_iterator pAend = corrframes[rID[0]].end();
#pragma omp parallel shared(cleanlist)// num_threads(20)
						{
#pragma omp for
	//for (Frame::const_iterator pA = corrframes[rID[0]].begin(); pA != pAend; ++pA) {
	for (int i = 0; i < corrframes[rID[0]].NumParticles(); i++) {
		Frame::const_iterator pA = corrframes[rID[0]].begin() + i;
//		if (fabs(pA->X() - 445.17794) < 1e-4) {
//			bool breakpoint = true;
//		}

		// this particle's coordinates in world space
		Position pAworld(cams[camID[0]].ImageToWorld(*pA));

		// 2nd camera parameters
		int Npixw = cams[camID[1]].Get_Npixw();
		int Npixh = cams[camID[1]].Get_Npixh();
		double kr = cams[camID[1]].Get_kr();
		double wpix = cams[camID[1]].Get_wpix();
		double hpix = cams[camID[1]].Get_hpix();
		double Noffh = cams[camID[1]].Get_Noffh();
		double Noffw = cams[camID[1]].Get_Noffw();
		//A structure that stores the cam parameters
		CamParam C1(Npixw, Npixh, kr, wpix, hpix, Noffh, Noffw);

		// position of 1st camera's projective center on second camera
		/* CHANGE 4/18/2013 */
		Position center0on1(cams[camID[1]].WorldToImage(cams[camID[0]].Center())); /* DIFF */
		Position particle0on1(cams[camID[1]].WorldToImage(pAworld));
		// unit vector in (projected) line of sight direction
		Position lineofsight(particle0on1 - center0on1);
		lineofsight /= lineofsight.Magnitude();
		// unit vector normal to the line of sight
		Position perpdir(lineofsight.Y(), -lineofsight.X(), 0); /* DIFF */

																// projected line: y = m*x + c (m = Line0on1.first, c = Line0on1.second)
		pair<double, double> Line0on1 = LineParam(camID[0], camID[1], pAworld);

		// A structure that stores the equation of lines at +-mindist_2D dist. from proj. line
		ParallelLines P_0on1(Line0on1.second, Line0on1.first, mindist_2D);

		// identifying particles on 2nd camera around the projected line (within mindist_2D) and pushing them into a tmplist
		vector<Frame::const_iterator> tmplist0on1;
		ParticleFinder1to1(C1, P_0on1, center0on1, perpdir, rID[1], tmplist0on1);

		avgsize += tmplist0on1.size();
		numlists++;

		if (rcams == 2) {
			throw runtime_error("Need at least 3 cameras for triangulation");
			continue;
		}


		// looping over each possible particle candidate in 2nd camera to find a match in subsequent cameras
		for (int it = 0; it < tmplist0on1.size(); it++) {
			Frame::const_iterator pB = tmplist0on1[it];
			Position pBworld(cams[camID[1]].ImageToWorld(*pB));
			vector<Frame::const_iterator> tmplist0n1on2;

			// list of particle candidates on 3rd camera close to the intersection of lines projected from 1st and 2nd cam
			ParticleFinder2to1(camID[2], rID[2], camID[0], camID[1], pAworld, pBworld, mindist_1D, tmplist0n1on2);

			pair<double, double> Line1on0 = LineParam(camID[1], camID[0], pBworld);

			// looping over each possible particle candidate in 3rd camera to find a match in subsequent cameras
			for (int itt = 0; itt < tmplist0n1on2.size(); itt++) {

				Frame::const_iterator pC = tmplist0n1on2[itt];
				Position pCworld(cams[camID[2]].ImageToWorld(*pC));

				// checking if the particle candidates in 2nd and 3rd cam are potential candidates
				// by projecting them back onto 1st and 2nd cams
				if (ParticleCheck2to1(camID[0], camID[1], camID[2], pBworld, pCworld, pA, mindist_1D) && ParticleCheck2to1(camID[1], camID[0], camID[2], pAworld, pCworld, pB, mindist_1D)) {
					if (rcams == 3) {
						try {
#pragma omp critical(cleanlist) // to indicate this sentence should only be executed by one thread at a time
							// to prevent conflict between different thread since it is shared.
							{
							cleanlist[0].push_back(pA); cleanlist[1].push_back(pB); cleanlist[2].push_back(pC);
							}
						}
						catch (bad_alloc) {
							runtime_error("cannot assign storage space for cleanlist");
						}
						//temps1++;
						continue;
					}

					// if ncams > 3, for 4th camera
					vector<Frame::const_iterator> tmplist1n2on3;
					// list of particle candidates on 4th camera close to the intersection of lines projected from 2nd and 3rd cam
					ParticleFinder2to1(camID[3], rID[3], camID[1], camID[2], pBworld, pCworld, mindist_1D, tmplist1n2on3);


					for (int ittt = 0; ittt < tmplist1n2on3.size(); ittt++) {

						Frame::const_iterator pD = tmplist1n2on3[ittt];
						// checking if the potential candidates on 4th cam are within mindist_2D to a line proj from 1st cam
						if (ParticleCheck1to1(camID[0], camID[3], pAworld, pD, mindist_1D)) {
							Position pDworld(cams[camID[3]].ImageToWorld(*pD));

							//checking if the particle candidates in 3rd and 4th cam are potential candidates
							// by projecting them back onto 1st and 2nd cams
							// edit: added all the combinations of checks to increase the robustness
							if (ParticleCheck2to1(camID[0], camID[2], camID[3], pCworld, pDworld, pA, mindist_1D) && ParticleCheck2to1(camID[0], camID[1], camID[3], pBworld, pDworld, pA, mindist_1D) && 
								ParticleCheck2to1(camID[1], camID[2], camID[3], pCworld, pDworld, pB, mindist_1D) && ParticleCheck2to1(camID[1], camID[0], camID[3], pAworld, pDworld, pB, mindist_1D) &&
								ParticleCheck2to1(camID[2], camID[0], camID[3], pAworld, pDworld, pC, mindist_1D) && ParticleCheck2to1(camID[2], camID[1], camID[3], pBworld, pDworld, pC, mindist_1D) &&
								ParticleCheck2to1(camID[3], camID[0], camID[1], pAworld, pBworld, pD, mindist_1D) && ParticleCheck2to1(camID[3], camID[0], camID[2], pAworld, pCworld, pD, mindist_1D)) {
								if (rcams >= 4) {
									try {
#pragma omp critical(cleanlist)
							{
										cleanlist[0].push_back(pA); cleanlist[1].push_back(pB); cleanlist[2].push_back(pC); cleanlist[3].push_back(pD);
							}
									}
									catch (bad_alloc) {
										runtime_error("cannot assign storage space for cleanlist");
									}

									if (rcams == 4)
										continue;

									// if ncams > 4, all the subsequent cams will only be used to check if the match from pA to pD is correct
									for (int Cam = 4; Cam < rcams; Cam++) {
										vector<Frame::const_iterator> tmplist2n3onCam;
										ParticleFinder2to1(camID[Cam], rID[Cam], camID[2], camID[3], pCworld, pDworld, mindist_1D, tmplist2n3onCam);
										// deleting the match if no particle is found
										if (tmplist2n3onCam.size() == 0) {
#pragma omp critical(cleanlist)
							{
											cleanlist[0].pop_back(); cleanlist[1].pop_back(); cleanlist[2].pop_back(); cleanlist[3].pop_back();
							}
										}
										else 
											for (int itttt = 0; itttt < tmplist2n3onCam.size(); itttt++) {
												Frame::const_iterator pCam = tmplist2n3onCam[itttt];
												if (ParticleCheck1to1(camID[0], camID[Cam], pAworld, pCam, mindist_1D) && ParticleCheck1to1(camID[1], camID[Cam], pBworld, pCam, mindist_1D)) 
													continue;
												
												else {
#pragma omp critical(cleanlist)
							{
													cleanlist[0].pop_back(); cleanlist[1].pop_back(); cleanlist[2].pop_back(); cleanlist[3].pop_back();

							}
												}
											}
										
									}
								}
							}
						}
					}
				}
			}
		}
	}
						}

	//########## TERSTING TRIANGULATION ###########
	//Load_cleanlist("C:/Users/aks5577/Google Drive/PHD/2017-Fall/Ash-STB/Release/", 1, corrframes_Pixel, corrframes, cleanlist);

	std::cout << "\tNumber of matches = " << cleanlist[0].size() << endl;
	//std::cout << "\tNumber of failed checks on cam0 = " << temps1 << endl;
	//std::cout << "\tMean pairlist size: " << static_cast<double>(avgsize)/static_cast<double>(numlists) << endl;
	
	deque<Position> matchedPos;
	deque< deque<int> > frame_indices;
	deque<double> raydists;
	//std::cout << "\tPerforming 3D triangulation..." << endl;
			// we can match it!
	unsigned long int cleanlist_num = cleanlist[0].size();
#pragma omp parallel shared(matchedPos, frame_indices, raydists, PosTouse) //num_threads(8)
						{
#pragma omp for
//	for (long int p = 0; p < cleanlist[0].size(); p++) {
	for (unsigned long int p = 0; p < cleanlist_num; p++) {
//		double percent = p/total_cleanlist;
		deque<Position> PosToMatch;
		deque<int> indices;
		for (int i = 0; i < rcams; ++i) {
			try {
				PosToMatch.push_back(*(cleanlist[i][p]));
			}
			catch (bad_alloc&) {
				throw runtime_error("Failure to allocate storage for PosToMatch");
			}
			indices.push_back(cleanlist[i][p].where());
		}
		pair<double,Position> wpos = WorldPosition(PosToMatch, ignoreCam);
#pragma omp critical
						{
		if (wpos.first < mindist_3D * mindist_3D && boundary_check.Check(wpos.second)) {
			// keep the matches which are inside the boundary and tolerant error
			try {
				matchedPos.push_back(wpos.second);
			}
			catch (bad_alloc&) {
				throw runtime_error("Failure to allocate storage for matchedPos");
			}
			
			frame_indices.push_back(indices);
			raydists.push_back(wpos.first);
			deque<Position> ttmp;
			for (int i = 0; i < rcams; ++i) {
				ttmp.push_back(*cleanlist[i][p]);
			}
			PosTouse.push_back(ttmp); // the sequences of frame_indices, raydists and PosTouse are corresponding to each other
		}
							}
	}
						}


	// output all the triangulated positions.
//	NumDataIO<double> pos_io;
//	pos_io.SetFilePath("/home/tanshiyong/Documents/Data/SinglePhase/SD0125/Tracks/pos.txt");
//	int num_pos = matchedPos.size();
//	pos_io.SetTotalNumber(num_pos * 4); // pos + tri-error
//	double* pos_e = new double[num_pos * 4];
//
//	for (int i = 0; i < num_pos; ++ i) {
//		pos_e[i * 4] = matchedPos[i].X();
//		pos_e[i * 4 + 1] = matchedPos[i].Y();
//		pos_e[i * 4 + 2] = matchedPos[i].Z();
//		pos_e[i * 4 + 3] = raydists[i];
//	}
//
//	pos_io.WriteData((double*) pos_e);
//	delete[] pos_e;

	cout << "\tNumber of matches within tri-error = " << PosTouse.size() << endl;



//	unsigned int n_match = frame_indices.size();
//	int* indices_vector = new int[n_match * 4];
//	for (unsigned int i = 0; i < n_match; ++i) {
//		indices_vector[i * 4] = frame_indices[i][0];
//		indices_vector[i * 4 + 1] = frame_indices[i][1];
//		indices_vector[i * 4 + 2] = frame_indices[i][2];
//		indices_vector[i * 4 + 3] = frame_indices[i][3];
//	}

//	NumDataIO<int> data_io;
//	data_io.SetTotalNumber(n_match * 4);
//	data_io.SetFilePath("/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/Triangulation/match_indices.txt");
//	data_io.WriteData(indices_vector);

///*
	unsigned int num_match =  matchedPos.size();
	unsigned int num_particle = 0;
//	double preference[num_match][2];
	double** preference = new double*[num_match];
	for (unsigned int i = 0; i < num_match; ++i) {
		preference[i] =  new double[2];
		preference[i][0] = 0; preference[i][1] = raydists[i]; //triangulation error
	}
	for (unsigned int i = 0; i < rcams; i++) {
		num_particle = corrframes[rID[i]].NumParticles();
		int* buffer = new int[num_particle]; // the buffer is used to save the match index with smaller error
		for (unsigned int j = 0; j < num_particle; j++) buffer[j] = -1; // initialize buffer with -1
		GroupAndPickMin(raydists, frame_indices, buffer, num_particle, num_match, i);
		for (unsigned int j = 0; j < num_particle; ++j) {
			if (buffer[j] == -1) continue;
//			preference[buffer[j]][1] = preference[buffer[j]][1] * preference[buffer[j]][0] + raydists[buffer[j]]
//				                       / (preference[buffer[j]][0] + 1);  // Why? for each match it has only one triangulation error. TODO
			preference[buffer[j]][0] = preference[buffer[j]][0] + 1;
		}
		delete[] buffer;
	}

	// Sorting the matches sequence according to its preference
	int* match_sequence = new int[num_match];
	for (int i = 0; i < num_match; ++i) match_sequence[i] = i;
	MergeSort(match_sequence, preference, num_match); //TODO: confirm it is ranging correctly


//	double error[num_match];
//	int frameindex[num_match][4];
//	int matchindex[num_match];
//	int preference_num[num_match];
//	for (int i = 0; i < num_match; ++i) {
//		int index = match_sequence[i];
//		error[i] = raydists[index];
//		frameindex[i][0] = frame_indices[index][0];frameindex[i][1] = frame_indices[index][1];
//		frameindex[i][2] = frame_indices[index][2];frameindex[i][3] = frame_indices[index][3];
//		matchindex[i] = index;
//		preference_num[i] = preference[index][0];
//	}
//
//			NumDataIO<double> error_io;
//			error_io.SetFilePath("/storage/home/sut210/work/Experiment/EXP7/error.txt");
//			error_io.SetTotalNumber(num_match);
//			error_io.WriteData((double*) error);
//
//			NumDataIO<int> index_io;
//			index_io.SetFilePath("/storage/home/sut210/work/Experiment/EXP7/frameindex.txt");
//			index_io.SetTotalNumber(num_match * 4);
//			index_io.WriteData((int*) frameindex);
//			index_io.SetFilePath("/storage/home/sut210/work/Experiment/EXP7/matchindex.txt");
//			index_io.SetTotalNumber(num_match);
//			index_io.WriteData((int*) matchindex);
//			index_io.SetFilePath("/storage/home/sut210/work/Experiment/EXP7/preferencenum.txt");
//						index_io.SetTotalNumber(num_match);
//						index_io.WriteData((int*) preference_num);


	int* fill_in[rcams];
	for (unsigned int i = 0; i < rcams; ++i) {
		num_particle = corrframes[rID[i]].NumParticles();
		fill_in[i] = new int[num_particle];
		for (unsigned int j = 0; j < num_particle; ++j) fill_in[i][j] = -1;
//		for (unsigned int j = 0; j < num_particle; ++j) fill_in[i][j] = 0;
	}
	bool* selection = new bool[num_match]; // to avoid segmentation fault
	for (unsigned int i = 0; i < num_match; ++i) {
		selection[i] =false;
	}
//	int iteration_times = 4; //total iteration times
//	for (int iteration = 0; iteration < iteration_times; ++ iteration) {
		for (unsigned int i = 0; i < num_match; ++i) {
			int sequence_math_index = match_sequence[i];  // get the new match index
			//selection[sequence_math_index] =false;
//			if (selection[sequence_math_index]) continue; // if one match has been used in the previous iteraction, then skip this one.
			int empty_times = 0;
			double sum_error = 0;
			int equal_times = 0;
			for (unsigned int j = 0; j < rcams; ++j) {  // to check each of the particles in the new match
				int particle_index = frame_indices[sequence_math_index][j];
				int match_index = fill_in[j][particle_index];

				// when there is no match for the checked particle
				if (match_index == -1) {
					empty_times = empty_times + 1;
					if (empty_times < rcams) {
						goto final_out;
					} else {  // all of the particles have no match
						for (unsigned int k = 0; k < rcams; ++k) {
							fill_in[k][frame_indices[sequence_math_index][k]] = sequence_math_index; // fill in the match index if it is all empty
						}
						selection[sequence_math_index] = true;
						break;
					}
				}

				// when there has already been a match for the checked particle
				// the preference # is less than ONE of the preference # of those filled match, then discard the new match
				if (preference[sequence_math_index][0] < preference[match_index][0]) {
					break;
				} else if (preference[sequence_math_index][0] == preference[match_index][0]) {
					// when the preference # is the same as the current preference # of those filled matches,
					// then we should calculate the average of the preference error of those filled matches
					equal_times = equal_times + 1;
					sum_error = sum_error + preference[match_index][1];
				}

				final_out:
				if (j < rcams - 1) continue; // continue to next camera when it haven't reached out to the final camera

				if (equal_times > 0) { // when there is no larger preference number,
					//if the preference error is larger than the averaged of those filled match, then discard the new match
					if (preference[sequence_math_index][1] > sum_error / equal_times) break;
				}

				// if the new match run through all of the checks above:
				// the preference # is at least equal or larger than the current preference #
				// the preference error is less than the current average preference error
				for (unsigned int k = 0; k < rcams; ++k) {
					match_index = fill_in[k][frame_indices[sequence_math_index][k]];
					if (match_index != -1) {
						selection[match_index] = false;
						for (unsigned int l = 0; l < rcams; ++l) {
							fill_in[k][frame_indices[match_index][k]] = -1;
						}
					}
					fill_in[k][frame_indices[sequence_math_index][k]] = sequence_math_index; // replace the match index
				}
				selection[sequence_math_index] = true;
			}
		}
//	}

	// second round: fill in 2D centers which haven't been filled
//	for (unsigned int i = 0; i < num_match; ++i) {
//		int sequence_math_index = match_sequence[i];  // get the new match index
//		if (!selection[sequence_math_index]) {
//			for (unsigned int j = 0; j < rcams; ++j) {
//				if (fill_in[j][frame_indices[sequence_math_index][j]] == -1) {
//					selection[sequence_math_index] = true;
//					for (unsigned int k = 0; k < rcams; ++k) {
//						fill_in[k][frame_indices[sequence_math_index][k]] = sequence_math_index;
//					}
//					break;
//				}
//			}
//		}
//	}

//	int max_select_num = 3; // the number of times of selection for each 2D center
//	double max_triangulatin_err = 0;
//	double averag_triangulation_err = 0;
//	// calculate the average triangulation error
//	for (unsigned int i = 0; i < num_match; ++i) {
//		averag_triangulation_err = averag_triangulation_err + preference[i][1];
//	}
//	averag_triangulation_err = averag_triangulation_err / num_match;
//	// calculate the standard deviation error for the triangulation error
//	double std_triangulation_err = 0;
//	for (unsigned int i = 0; i < num_match; ++i) {
//		std_triangulation_err = std_triangulation_err + pow((preference[i][1] - averag_triangulation_err), 2);
//	}
//	std_triangulation_err = sqrt(std_triangulation_err / num_match);
//
//	max_triangulatin_err = averag_triangulation_err + 4 * std_triangulation_err;
//
//	for (unsigned int i = 0; i < num_match; ++i) {
//		int sequence_match_index = match_sequence[i];  // get the new match index
//		selection[sequence_match_index] =false;
//		bool availability = true;
//		for (unsigned int j = 0; j < rcams; ++j) {  // to check each of the particles in the new match
//			int particle_index = frame_indices[sequence_match_index][j];
//			if (fill_in[j][particle_index] > max_select_num) {
//				availability =  false; // indicate no more available 2D centers for this match
//			}
//		}
//		if (availability) {
//			if (preference[sequence_match_index][1] < max_triangulatin_err) {
//				selection[sequence_match_index] =true;
//				for (unsigned int j = 0; j < rcams; ++j) {  // to check each of the particles in the new match
//					int particle_index = frame_indices[sequence_match_index][j];
//					fill_in[j][particle_index] = fill_in[j][particle_index] + 1;
//				}
//			}
//		}
//	}


	for (unsigned int i = 0; i < num_match; ++i) delete[] preference[i];
	for (int i = 0; i < rcams; ++i) delete[] fill_in[i];
	delete[] match_sequence;
	delete[] preference;

//	vector<int> goodindex;

//	unsigned int num_match =  matchedPos.size();
	deque<Position> goodPos;
	int match_index = 0;
	for (unsigned int i = 0; i < num_match; ++i) {
		if (!selection[i]) continue;

//		goodindex.push_back(i);

		deque< deque<double> > pos2D(4);
		for (int id = 0; id < 4; id++) {
			if (id < rcams) {
				deque<double> tmp2D(2);
//				Position temp = cams[id].Distort((PosTouse[i])[id]);
				Position temp = cams[camID[id]].Distort((PosTouse[i])[id]);
				tmp2D[0] = temp.X(); tmp2D[1] = temp.Y();
				pos2D[id] = tmp2D;
			} else {
				deque<double> tmp(2);
				tmp[0] = 0.0; tmp[1] = 0.0;
				pos2D[id] = tmp;
			}
		}
		Position worldposi(matchedPos[i].X(), matchedPos[i].Y(), matchedPos[i].Z(),
							pos2D[0][0], pos2D[0][1], pos2D[1][0], pos2D[1][1], pos2D[2][0], pos2D[2][1], pos2D[3][0], pos2D[3][1], 0);
		good2Dpos.push_back(worldposi);
		goodPos.push_back(matchedPos[i]);
	}

	delete[] selection;
//	*/
//	pos_io.SetFilePath("/home/tanshiyong/Documents/Data/SinglePhase/SD0125/Tracks/select_pos.txt");
//	num_pos = goodPos.size();
//	pos_io.SetTotalNumber(num_pos * 3); // pos + tri-error
//	double* pos_2 = new double[num_pos * 3];
//
//	for (int i = 0; i < num_pos; ++ i) {
//		pos_2[i * 3] = goodPos[i].X();
//		pos_2[i * 3 + 1] = goodPos[i].Y();
//		pos_2[i * 3 + 2] = goodPos[i].Z();
//	}
//
//	pos_io.WriteData((double*) pos_2);
//	delete[] pos_2;

/*
	deque<Position> goodPos;
	for (unsigned int i = 0; i < PosTouse.size(); ++i) {
//		if (!selection[i]) continue;

//		goodindex.push_back(i);

		deque< deque<double> > pos2D(4);
		for (int id = 0; id < 4; id++) {
			if (id < rcams) {
				deque<double> tmp2D(2);
//				Position temp = cams[id].Distort((PosTouse[i])[id]);
				Position temp = cams[camID[id]].Distort((PosTouse[i])[id]);
				tmp2D[0] = temp.X(); tmp2D[1] = temp.Y();
				pos2D[id] = tmp2D;
			} else {
				deque<double> tmp(2);
				tmp[0] = 0.0; tmp[1] = 0.0;
				pos2D[id] = tmp;
			}
		}
		Position worldposi(matchedPos[i].X(), matchedPos[i].Y(), matchedPos[i].Z(),
							pos2D[0][0], pos2D[0][1], pos2D[1][0], pos2D[1][1], pos2D[2][0], pos2D[2][1], pos2D[3][0], pos2D[3][1], 0);
		good2Dpos.push_back(worldposi);
		goodPos.push_back(matchedPos[i]);
	}
	*/
//	int num_good = goodindex.size();
//		double error[num_good];
//		int frameindex[num_good][4];
//		int matchindex[num_good];
//		for (int i = 0; i < num_good; ++i) {
//			int index = goodindex[i];
//			error[i] = raydists[index];
//			frameindex[i][0] = frame_indices[index][0];frameindex[i][1] = frame_indices[index][1];
//			frameindex[i][2] = frame_indices[index][2];frameindex[i][3] = frame_indices[index][3];
//			matchindex[i] = index;
//		}
//
//
//		NumDataIO<double> error_io;
//		error_io.SetFilePath("/storage/home/sut210/work/Experiment/EXP7/error.txt");
//		error_io.SetTotalNumber(num_good);
//		error_io.WriteData((double*) error);
//
//		NumDataIO<int> index_io;
//		index_io.SetFilePath("/storage/home/sut210/work/Experiment/EXP7/frameindex.txt");
//		index_io.SetTotalNumber(num_good * 4);
//		index_io.WriteData((int*) frameindex);
//		index_io.SetFilePath("/home/sut210/Documents/Experiment/EXP2/matchindex.txt");
//		index_io.SetTotalNumber(num_good);
//		index_io.WriteData((int*) matchindex);

//	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//cout << "\tmindist_2D: " << mindist_2D << endl;
	cout << "\tNumber of triangulated particles = " << goodPos.size() << endl;
//	std::cout << "\tMatching time: " << duration << endl;

	delete[] cleanlist;

	return Frame(goodPos);
}

int Calibration::GroupAndPickMin( deque<double>& raydists, deque< deque<int> >& frame_index, int* buffer,
		int num_particle, int num_match, int camera_num) {

	int particle_index = 0;
//	if (is_list_empty) {
		for (int i = 0; i < num_match; ++i) {
			particle_index = frame_index[i][camera_num];
			if (buffer[particle_index] == -1) {
				buffer[particle_index] = i;
			} else {
				if (raydists[buffer[particle_index]] > raydists[i]) {
					buffer[particle_index] = i; // replace the buffer with a match with smaller error
				}
			}
		}
//	} else {
//		int match_index = 0;
//		for (int i = 0; i < list_size; ++i) {
//			match_index = minimum_list[i];
//			if (match_index == -1) continue; // skip those empty index;
//			particle_index = frame_index[match_index][camera_num]; // get the particle index of the match in the specific camera
//			if (buffer[particle_index] == -1) {
//					buffer[particle_index] = match_index;
//			} else {
//					if (raydists[buffer[particle_index]] > raydists[match_index]) {
//						buffer[particle_index] = match_index; // replace the buffer with a match with smaller error
//					}
//			}
//		}
//	}

	return 0;
}

pair<double,Position> Calibration::WorldPosition(deque<Position> ipos, int ignoreCam) throw(runtime_error)
{
	/*if (ipos.size() != cams.size()) {
		throw runtime_error("Number of cameras and number of images do not match!");o
	}*/

	// solve, in a least-squares sense, the intersection of the three lines of sight.
	// this will give us the distance of closest approach to the three lines.
	// embarrassingly, see
	// http://en.wikipedia.org/wiki/Line-line_intersection

	Matrix M;
	Position P(0, 0, 0);

	//Position sight[ncams];
    Position sight[20];

	for (int i = 0; i < rcams; ++i) {
		// construct a line of sight for this point:
		// vector from camera center to this point (in world coordinates)
		sight[i] = cams[camID[i]].ImageToWorld(ipos[i]) - cams[camID[i]].Center();
		// normalize
		sight[i] /= sight[i].Magnitude();

		// add to the least-squares matrices

		// inelegant, but it works
		Matrix tmp;
		tmp.Set(0, 0, 1 - sight[i].X() * sight[i].X());
		tmp.Set(0, 1, -sight[i].X() * sight[i].Y());
		tmp.Set(0, 2, -sight[i].X() * sight[i].Z());
		tmp.Set(1, 0, -sight[i].Y() * sight[i].X());
		tmp.Set(1, 1, 1 - sight[i].Y() * sight[i].Y());
		tmp.Set(1, 2, -sight[i].Y() * sight[i].Z());
		tmp.Set(2, 0, -sight[i].Z() * sight[i].X());
		tmp.Set(2, 1, -sight[i].Z() * sight[i].Y());
		tmp.Set(2, 2, 1 - sight[i].Z() * sight[i].Z());

		P += tmp * cams[camID[i]].Center();
		M += tmp;
	}

	// invert the matrix and construct the 3D position
	Position tmpi(M.Invert() * P);
	Position worldpos(M.Invert() * P);

	// calculate the rms distance from worldpos to each ray
	double dist = 0;
	for (int i = 0; i < rcams; ++i) {
		Position h = (worldpos - Dot(worldpos, sight[i]) * sight[i]
									- (cams[camID[i]].Center() - Dot(cams[camID[i]].Center(), sight[i]) * sight[i]));
		dist += h.Magnitude2();
	}
	dist /= static_cast<double>(rcams);

	deque< deque<double> > pos2D(4);
	for (int i = 0; i < 4; i++) {
		if (i < ncams && i != ignoreCam) {
			deque<double> tmp2D(2);
			Position tmp = cams[i].WorldToImage(tmpi); tmp.Set_Z(0);
			tmp = cams[i].Distort(tmp);
			tmp2D[0] = tmp.X(); tmp2D[1] = tmp.Y();
			pos2D[i] = tmp2D;
		}
		else {
		/*
		 * Modified by Shiyong Tan, 1/30/18
		 * Initialization tmp = {0.0, 0.0} is not permited by c++98
		 * Use tmp(0.0, 0.0) instead
		 * Start;
		 */
//				deque<double> tmp = { 0.0, 0.0 };
			deque<double> tmp(2); //TODO: Check such initialization is valid.
			tmp[0] = 0.0; tmp[1] = 0.0;
		// End
			pos2D[i] = tmp;
		}
	}
	Position worldposi(tmpi.X(), tmpi.Y(), tmpi.Z(), pos2D[0][0], pos2D[0][1], pos2D[1][0], pos2D[1][1], pos2D[2][0], pos2D[2][1], pos2D[3][0], pos2D[3][1], dist);
	return make_pair(dist, worldposi);	
}

// a function that returns the slope (m_Proj) and intercept (c_Proj) of the line projected by camid2 on camid1
pair<double, double> Calibration::LineParam(int camid1, int camid2, Position P) 
{
	Position center(cams[camid2].WorldToImage(cams[camid1].Center()));
	Position particle(cams[camid2].WorldToImage(P));
	double m_Proj = (center.Y() - particle.Y()) / (center.X() - particle.X());
	double c_Proj = particle.Y() - m_Proj*particle.X();
	return make_pair(m_Proj, c_Proj);
}

// A function thats takes all particle positions as an input and gives the pixel matrix is_Particle (with 1 when a particle is present and 0 when there are no particles). 
// It also gives is_ParticlePos (a 2D array of deuqe<Position>) that stores all particle centers when is_Particle is 1. 
void Calibration::Peak_identifier(int camID, int rID, Frame& f, Frame& f2)
{
		int Npixw = cams[camID].Get_Npixw();
		int Npixh = cams[camID].Get_Npixh();

		// setting is_Peak elements to false and is_ParticlePos deques to empty
		for (int i = 0; i < Npixh; i++)	{
			memset(Calibration::is_Particle[rID][i], false, Npixh * sizeof(bool));
			for (int j = 0; j < Npixw; j++)	{
				Calibration::is_ParticlePos[rID][i][j].clear();
			}
		}

		Frame::const_iterator pend = f.end();
		Frame::const_iterator p2end = f2.end();
		Frame::const_iterator t = f.begin();
		for (Frame::const_iterator t2 = f2.begin(); t2 != p2end; ++t2)	{
			int k = (*t).X();
			int l = (*t).Y();

			if (k > Npixw-1 || l > Npixh-1 || k < 0 || l < 0)	{
				//cout << "particle out of image bounds: " << l << "," << k << endl;
				++t;
				continue;
			}
			Calibration::is_Particle[rID][l][k] = true;
			Calibration::is_ParticlePos[rID][l][k].push_back(t2);
			++t;
	}
}

/* Steps to find the particles within mindist_2D distance from the projected line,
Step 1: Identify 2 parallel lines at +-midist_2D distance from projected line. (struct ParallelLines)
Step 2: Obtain the intersection of the 2 lines with image pixels_orig. (Functions PixelSearchx or PixelSearchy)
		This creates a search bandwidth of pixels_orig around the projected line. 
Step 3: Use is_Particle logic matrix to check if the pixels_orig in this badwidth have a particle, 
		and then push the particles into pairlist. (Function ParticleFinder)
*/

// Step 1:
// The structure ParallelLines defined in the begining of this class takes care of this step.

// Step 2:
// PixelSearchx or PixelSearchy to identify the search bandwidth of pixels around the projected line (y = mx + c in physical units(mm)).
// PixelSearchx and PixelSearchy are used when abs(slope,m) < 1 and abs(slope,m) > 1 respectively.

/* Note: Due to distortion, the projected straight line in mm appears as a curved line in pixels. 
		 Therefore, the projected line (y=mx+c) is always in physical units (mm) */
void Calibration::PixelSearchx(CamParam C, ParallelLines P, double ypix, int& Min, int& Max)
{
	int temp = 0;	
	double x_a, x_b, x_c;

	if (P.m_proj > 0)	{
		x_a = Xpix1(C, P.m_proj, P.c_proj, ypix);
		x_b = Xpix1(C, P.m_proj, P.c_plusdel, ypix);
		x_c = Xpix1(C, P.m_proj, P.c_minusdel, ypix);
	}

	if (P.m_proj < 0)	{
		x_a = Xpix2(C, P.m_proj, P.c_proj, ypix);
		x_b = Xpix2(C, P.m_proj, P.c_plusdel, ypix);
		x_c = Xpix2(C, P.m_proj, P.c_minusdel, ypix);
	}

	int XPix_a = x_a;
			// only considering pixels inside the image bounds
	if (XPix_a >= 0 || XPix_a < C.Npixw)	{
		temp = x_b;
		temp = max(1, temp);
		int YPix_b = min(temp, C.Npixw - 1);
		temp = x_c;
		temp = max(1, temp);
		int YPix_c = min(temp, C.Npixw - 1);

		Min = min(YPix_b, YPix_c)-1;
		Max = max(YPix_b, YPix_c);
		Max = min(Max, C.Npixw - 2)+2;
	}

}

void Calibration::PixelSearchy(CamParam C, ParallelLines P, double xpix, int& Min, int& Max)

{
	int temp = 0;
	double y_a, y_b, y_c;

	y_a = Ypix1(C, P.m_proj, P.c_proj, xpix);
	y_b = Ypix1(C, P.m_proj, P.c_plusdel, xpix);
	y_c = Ypix1(C, P.m_proj, P.c_minusdel, xpix);
	
	int YPix_a = y_a;
		// only considering the pixels_orig inside the image bounds
	if (YPix_a >= 0 || YPix_a < C.Npixh)	{
		temp = y_b;
		temp = max(1, temp);
		int YPix_b = min(temp, C.Npixh - 1);
		temp = y_c;
		temp = max(1, temp);
		int YPix_c = min(temp, C.Npixh - 1);

		Min = min(YPix_b, YPix_c) - 1;
		Max = max(YPix_b, YPix_c);
		Max = min(Max, C.Npixh - 2) + 2;
	}
}

/* 
Due to radial distortion, the intersection of lines y = mx + c (mm) with image pixels results in a quadratic equation.
=>Functions Ypix1 and Ypix2 return the roots of the quadratic equation when scanning in x-direction.
=>Functions Xpix1 and Xpix2 return the roots of the quadratic equation when scanning in y-direction.
*/
double Calibration::Ypix1(CamParam P, double m, double c, double Xpix)
{
	double X = (Xpix - P.Npixw/2 + P.Noffw)*P.wpix;
	double a, y;
	if (P.kr == 0)
		y = m*X + c;
	else {
		a = c*sqrt(1 - 4 * pow(X, 2) * P.kr*(1 + pow(m, 2) + P.kr*pow(c, 2)) - 4 * X*c*P.kr*m);
		y = (c + 2 * X * m + a) / (2 * (P.kr * pow(c, 2) + 1));
	}
	y = (P.Npixh / 2 - P.Noffh) - (y / P.hpix);
	return y;
}

double Calibration::Ypix2(CamParam P, double m, double c, double Xpix)
{
	double X = (Xpix - P.Npixw / 2 + P.Noffw)*P.wpix;
	double a, y;
	if (P.kr == 0)
		y = m*X + c;
	else {
		a = c*sqrt(1 - 4 * pow(X, 2) * P.kr*(1 + pow(m, 2) + P.kr*pow(c, 2)) - 4 * X*c*P.kr*m);
		y = (c + 2 * X * m - a) / (2 * (P.kr * pow(c, 2) + 1));
	}
	y = (P.Npixh / 2 - P.Noffh) - (y / P.hpix);
	return y;
}

double Calibration::Xpix1(CamParam P, double m, double c, double Ypix)
{
	double Y = (-Ypix + P.Npixh / 2 - P.Noffh)*P.hpix;
	double a, x;
	if (P.kr == 0)
		x = -c + Y / m;
	else {
		a = c*sqrt(pow(m, 2) - 4 * pow(Y, 2) * P.kr*(1 + pow(m, 2) + P.kr*pow(c, 2)) + 4 * Y*c*P.kr*m);
		x = (-c*m + 2 * Y * m - a) / (2 * (P.kr * pow(c, 2) + pow(m, 2)));
	}
	x = (P.Npixw / 2 - P.Noffw) + (x / P.wpix);
	return x;
}

double Calibration::Xpix2(CamParam P, double m, double c, double Ypix)
{
	double Y = (-Ypix + P.Npixh / 2 - P.Noffh)*P.hpix;
	double a, x;
	if (P.kr == 0)
		x = -c + Y / m;
	else {
		a = c*sqrt(pow(m, 2) - 4 * pow(Y, 2) * P.kr*(1 + pow(m, 2) + P.kr*pow(c, 2)) + 4 * Y*c*P.kr*m);
		x = (-c*m + 2 * Y * m + a) / (2 * (P.kr * pow(c, 2) + pow(m, 2)));
	}
	x = (P.Npixw / 2 - P.Noffw) + (x / P.wpix);
	return x;
}


// Step 3: 
// ParticleFinder: A function that finds the particles in the pixel search bandwidth and pushes them into tmplist.
void Calibration::ParticleFinder1to1(CamParam C,ParallelLines P, Position center, Position perpdir, int k, vector<Frame::const_iterator>& tmplist)
{
	int Npix = (abs(P.m_proj) >= 1) ? C.Npixh : C.Npixw;
	// scanning through the image pixels in x or y direction depending on the slope of the line
	for (int pix = 0; pix < Npix; pix++)
	{
		int xPix, yPix;
		int Min = 0;
		int Max = 0;

		// if slope of proj line is greater than 1, scanning in y - direction
		if (abs(P.m_proj) >= 1)
		{
			yPix = pix;
			PixelSearchx(C, P, pix, Min, Max);
			for (xPix = Min; xPix < Max; xPix++)
			{
				// checking for pixel candidates with particles
				if (is_Particle[k][yPix][xPix] == 1) {
					// storing the particles in that pixel
					for (int it = 0; it < is_ParticlePos[k][yPix][xPix].size(); it++) {
						Frame::const_iterator pB = is_ParticlePos[k][yPix][xPix][it];
						// if the distance from camera i's projective center along perpdir
						// is less than mindist_2D, this is a potential match!
						Position pBline(*pB - center);
						if (abs(Dot(pBline, perpdir)) < P.mindist_2D)
							tmplist.push_back(pB);
					}
				}
			}
		}

		// if slope of proj line is less than 1, scanning in x - direction
		else if (abs(P.m_proj) < 1)	{
			PixelSearchy(C, P, pix, Min, Max);
			xPix = pix;
			for (yPix = Min; yPix < Max; yPix++) {
				if (is_Particle[k][yPix][xPix] == 1) {
					for (int it = 0; it < is_ParticlePos[k][yPix][xPix].size(); it++) {
						Frame::const_iterator pB = is_ParticlePos[k][yPix][xPix][it];
						Position pBline(*pB - center);
						if (abs(Dot(pBline, perpdir)) < mindist_2D)
							tmplist.push_back(pB);
					}
				}
			}
		}
	}
}

// A function that creates a search area around the intersection of 2 projected lines (from camID 1 and camid2) on camID
// and pushes all the particles within this search area to tmplist
void Calibration::ParticleFinder2to1(int camID, int rID, int camid1, int camid2, Position P1world, Position P2world, double mindist_1D, vector<Frame::const_iterator>& tmplist)
{
	// projecting the lines from camid1, camid2 onto camID and finding their intersection (x,y)
	Position center1(cams[camID].WorldToImage(cams[camid1].Center()));
	Position particle1(cams[camID].WorldToImage(P1world));
	double m1 = (center1.Y() - particle1.Y()) / (center1.X() - particle1.X());
	double c1 = particle1.Y() - m1*particle1.X();
	Position center2(cams[camID].WorldToImage(cams[camid2].Center()));
	Position particle2(cams[camID].WorldToImage(P2world));
	double m2 = (center2.Y() - particle2.Y()) / (center2.X() - particle2.X());
	double c2 = particle2.Y() - m2*particle2.X();
	double x = (c2 - c1) / (m1 - m2);
	double y = m1*x + c1;

	// getting the acute angle between the 2 lines
	double theta1 = atan(m1)*(180 / PI), theta2 = theta2 = atan(m2) * (180.0 / PI);
	double angle = abs(theta1 - theta2);
	if (angle > 180)
		angle = 360 - angle;

//	double mindist = abs(mindist_1D / sin(PI*angle / 360));
//	if (mindist > 1.5)
//		mindist = 1.5;
	double mindist = abs(mindist_1D / sin(PI*angle / 360));
	/*
	 * Modified by Shiyong Tan
	 * Date: 6.7.18
	 * Increase the threshold to eliminate the gap
	 * Start:
	 */
	if (mindist > 1200 * config.factor) //1200 voxels
		mindist = 1200 * config.factor;
	// End

	int Npixw = cams[camID].Get_Npixw();
	int Npixh = cams[camID].Get_Npixh();
	Position p(x, y, 0);
	// forming a search area around the intersection pointExtra
	Position p1(x - mindist, y - mindist, 0);
	Position p2(x + mindist, y + mindist, 0);
	Position pPix1(cams[camID].Distort(p1));
	Position pPix2(cams[camID].Distort(p2));
	// bounds of search area in pixels_orig
	int xpix1 = floor(pPix1.X());
	int xpix2 = ceil(pPix2.X());
	int ypix1 = ceil(pPix1.Y());
	int ypix2 = floor(pPix2.Y());

	// if any of those pixels are within the image bounds, check for particles in those pixels_orig
	if (((xpix1 >= 0 && xpix1 < Npixw) || (xpix2 >= 0 && xpix2 < Npixw)) || ((ypix1 >= 0 && ypix1 < Npixh) || (ypix2 >= 0 && ypix2 < Npixh)) || (xpix1*xpix2<=0 && ypix1*ypix2<=0))
	{
		xpix1 = max(0, xpix1);
		ypix2 = max(0, ypix2);
		xpix2 = min(Npixw-1, xpix2);
		ypix1 = min(Npixh-1, ypix1);
		for (int xpix = xpix1; xpix <= xpix2; xpix++)	
			for (int ypix = ypix2; ypix <= ypix1; ypix++)	
				if (is_Particle[rID][ypix][xpix] == 1)	
					for (int it = 0; it < is_ParticlePos[rID][ypix][xpix].size(); it++)	{
						Frame::const_iterator pX = is_ParticlePos[rID][ypix][xpix][it];			
						Position dist(*pX - p);
						if (dist.Magnitude() < mindist) {						// if the distance of the particle from intersection point is smaller than midist and
							double dist1 = abs(-m1*pX->X() + pX->Y() - c1) / sqrt(1 + pow(m1, 2));	// distance of particle from line1 
							double dist2 = abs(-m2*pX->X() + pX->Y() - c2) / sqrt(1 + pow(m2, 2));	// distance of particle from line2
							if (dist1 < mindist_1D && dist2 < mindist_1D) {	// if the particle is within mindist_1D dist from both the lines,
								//cout << mindist_1D << "," << mindist << "," << angle << "," << sin(PI*angle/360) << endl;
								tmplist.push_back(pX);							// add it to the list
							}
						}
					}	
	}
}

// A function to check if a particle pair (P1, P2) on camid1 and camid2 are a potential match for particle P0 on camid0
bool Calibration::ParticleCheck2to1(int camid0, int camid1, int camid2, Position P1world, Position P2world, Frame::const_iterator P0, double mindist_1D)
{
	// projecting the lines from camid1, camid2 onto camID and finding their intersection (x,y)
	Position center1(cams[camid0].WorldToImage(cams[camid1].Center()));
	Position particle1(cams[camid0].WorldToImage(P1world));
	double m1 = (center1.Y() - particle1.Y()) / (center1.X() - particle1.X());
	double c1 = particle1.Y() - m1*particle1.X();

	Position center2(cams[camid0].WorldToImage(cams[camid2].Center()));
	Position particle2(cams[camid0].WorldToImage(P2world));
	double m2 = (center2.Y() - particle2.Y()) / (center2.X() - particle2.X());
	double c2 = particle2.Y() - m2*particle2.X();
	double x = (c2 - c1) / (m1 - m2);
	double y = m1*x + c1;

	// getting the acute angle between the 2 lines
	double theta1 = atan(m1)*(180 / PI), theta2 = theta2 = atan(m2) * (180.0 / PI);
	double angle = abs(theta1 - theta2);
	if (angle > 180)
		angle = 360 - angle;

	/*
	 * Modified by Shiyong Tan
	 * Date: 6.7.18
	 * Increase the threshold to eliminate the gap
	 * Start:
	 */
	double mindist = abs(mindist_1D / sin(PI*angle / 360));
	if (mindist > 1200 * config.factor) //1200 voxels
		mindist = 1200 * config.factor;
	// END

	bool isMatch = false;
	double dist = pow((x - (*P0).X()), 2) + pow((y - (*P0).Y()), 2);
	if (dist <= pow(mindist, 2)) {
		double dist1 = abs(-m1*P0->X() + P0->Y() - c1) / sqrt(1 + pow(m1, 2));	// distance of particle from line1 
		double dist2 = abs(-m2*P0->X() + P0->Y() - c2) / sqrt(1 + pow(m2, 2));	// distance of particle from line2
		if (dist1 < mindist_1D && dist2 < mindist_1D)
			isMatch = true;
	}
	return isMatch;
}

// A function to check if a particle P1 on camid1 is a potential match for particle P0 on camid0
bool Calibration::ParticleCheck1to1(int camid0, int camid1, Position P0world, Frame::const_iterator P1, double mindist_1D) {
	bool isMatch = false;
	Position center(cams[camid1].WorldToImage(cams[camid0].Center()));
	Position particle(cams[camid1].WorldToImage(P0world));
	Position line(particle - center);
	line /= line.Magnitude();
	Position perp(line.Y(), -line.X(), 0);
	Position p1line(*P1 - center);

	if (abs(Dot(p1line, perp)) < mindist_1D)
		isMatch = true;

	return isMatch;
}

void Calibration::Load_cleanlist(string path, int frame, deque<Frame>& corrFrames_pixel, deque<Frame>& corrFrames, vector<Frame::const_iterator>*& cleanlist) {

/*
 * Modified by Shiyong Tan, 2/6/18
 * Discard using matio, use DataIO instead.
 * Start:
 */
	
//	cout << "break1" << endl;
//	stringstream s; s << path << "random_updatedCams_Spherical100000" << ".mat";
//	string file = s.str();
//	const char *fileName = file.c_str();
//	mat_t *mat = Mat_Open(fileName, MAT_ACC_RDONLY);
//
//	if (mat == NULL)
//		cout << " error in reading the cleanlist mat file" << endl;
//
//	stringstream ss; ss << "frame" << 1;
//	string varName = ss.str();
//
//	if (mat) {
//
//		matvar_t *matVar = 0;
//		matVar = Mat_VarRead(mat, (char*)varName.c_str());
//
//		if (matVar) {
//			int rows;	int cols;
//			rows = matVar->dims[0]; cols = matVar->dims[1];
//			unsigned namesize = matVar->nbytes / matVar->data_size;
//			double *namedata = static_cast<double*>(matVar->data);
//			for (int id = 0; id < 4; id++) {
//				Frame pos2D, pos2D_pixel;
//				for (int i = 0; i < 10000; i++) {
//					Position pos(namedata[(2 * id + 3) * rows + i], namedata[(2 * id + 4) * rows + i], 0);
//					pos2D_pixel.Add(pos);
//					pos2D.Add(cams[id].UnDistort(pos));
//				}
//				corrFrames.push_back(pos2D);
//				corrFrames_pixel.push_back(pos2D_pixel);
//			}
//		}
//	}
//	else
//		cout << "Cannot open file\n";

	//TODO: check whether it works. Shiyong Tan
	NumDataIO<double> data_io;
	data_io.SetFilePath(path + "random_updatedCams_Spherical100000.txt"); // Read txt file.
	//data format: 3D position(X,Y,Z) + 2D position(X,Y) for 4 Cameras, 10000 points
	int num_camera = 4; int num_points = 10000;
	int cols = 3 + 2 * num_camera; int rows = num_points;
	data_io.SetTotalNumber(cols * rows);
	double points_array[rows][cols];
	data_io.ReadData((double*) points_array);

	// Put points_array into pos
	for (int id = 0; id < num_camera; id++) {
		Frame pos2D, pos2D_pixel;
		for (int i = 0; i < num_points; i++) {
			Position pos(points_array[2 * id + 3][i], points_array[2 * id + 4][i], 0);
			pos2D_pixel.Add(pos);
			pos2D.Add(cams[id].UnDistort(pos));
		}
		corrFrames.push_back(pos2D);
		corrFrames_pixel.push_back(pos2D_pixel);
	}
// End

	Frame::const_iterator pA = corrFrames[0].begin();
	Frame::const_iterator pB = corrFrames[1].begin();
	Frame::const_iterator pC = corrFrames[2].begin();
	Frame::const_iterator pD = corrFrames[3].begin();

	for (int i = 0; i < corrFrames[0].NumParticles(); i++) {
		cleanlist[0].push_back(pA); cleanlist[1].push_back(pB); cleanlist[2].push_back(pC); cleanlist[3].push_back(pD);
		++pA; ++pB; ++pC; ++pD;
	}
//
//
//	Mat_Close(mat);
}

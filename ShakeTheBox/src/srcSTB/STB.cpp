

#include <STB.h>
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <utility> // to use smart pointer
#include "Position.h"
#include "NumDataIO.h"
#include "Common.h"

#define ALL_CAMS 100
#define REDUCED_CAMS 0
#define STBflag true
#define IPRflag false

using namespace std;

//enum STB::TrackType { Inactive = 0, ActiveShort = 1, ActiveLong = 2, Exit = 3, InactiveLong = 4};

// Remove the 'temp' frame after testing

STB::STB(int firstFrame, int lastFrame, string pfieldfile, string iprfile, int ncams, deque<int> camIDs, 
	deque<string> imageNameFiles, double initialRadius, double avgDist, double expShift, double masc, 
	double mrsc, double fpt, double lowerInt, bool iprFlag) : 

	first(firstFrame), last(lastFrame), iprfile(iprfile), _ipr(iprfile, ncams), pfieldfile(pfieldfile), ncams(ncams), camID(camIDs),
	imgNameFiles(imageNameFiles), searchRadiusSTB(initialRadius), avgIPDist(avgDist), largestShift(expShift), maxAbsShiftChange(masc), 
	maxRelShiftChange(mrsc/100), iprFlag(iprFlag), OTFcalib(ncams, _ipr.otfFile), fpt(fpt), intensityLower(lowerInt)
{ 
// SETTING UP THE PARAMETERS
	// camera parameters
	cams = _ipr.camsAll;
	Npixh = cams[0].Get_Npixh();
	Npixw = cams[0].Get_Npixw();

	// creating & initializing original, residual and reprojected image pixels
	for (int n = 0; n < ncams; n++) {
		pixels_orig.push_back(new int*[Npixh]);
		pixels_reproj.push_back(new int*[Npixh]);
		pixels_res.push_back(new int*[Npixh]);

		for (int i = 0; i < Npixh; i++) {
			pixels_orig[n][i] = new int[Npixw] {};
			pixels_reproj[n][i] = new int[Npixw] {};
			pixels_res[n][i] = new int[Npixw] {};
		}
	}

//	for (int i = 0; i < 100; ++i) {
//		linear_fit_error[i] = 0;
//	}

	it_innerloop = _ipr.it_innerloop;
	tiffaddress = _ipr.tiffaddress;
	
	for (int i = 0; i < ncams; i++) {												// loading image filenames at each frame
		ifstream file(tiffaddress + imgNameFiles[i], ios::in);
		if (file.is_open()) {
			deque<string> imgNames;
			int lineNum = 0;
			for (int frame = 1; frame <= last; frame++) {
				string line;
				getline(file, line); //line.erase(line.find_first_of(' '));
				stringstream tiffname;
				tiffname << tiffaddress << line;
				imgNames.push_back(tiffname.str());
			}
			imgSequence.push_back(imgNames);
		}
		else
			cout << "could not open the image sequence file " << imageNameFiles[i] << endl;
	}

	string address = tiffaddress + "Tracks/InitialTracks/";
	if ( !((debug_mode == SKIP_PREVIOUS_TRACKS || debug_mode == SKIP_PREVIOUS_BACK_STB) && debug_frame_number >= 4)) {
		if (debug_mode == SKIP_INITIAL_PHASE ) {
			LoadAllTracks(address, to_string(firstFrame + 3), 0);
		} else {
			InitialPhase(pfieldfile);														// Initial phase: obtaining tracks for first 4 time steps
			MatTracksSave(address, to_string(firstFrame + 3), 0);
		}
	}



	clock_t start0;
	start0 = clock();
	if (!(_ipr.triangulationOnly) && !(_ipr.IPROnly)) {								// if the user chooses to run only triangulation or IPR, stop here

		// Loading tracks from mat files
	/*	TrackType types;
		types = static_cast<TrackType>(0);
		Load_Tracks("InactiveTracks", types);
		types = static_cast<TrackType>(1);
		Load_Tracks("ActiveShortTracks", types);
		types = static_cast<TrackType>(2);
		Load_Tracks("ActiveLongTracks", types);
		types = static_cast<TrackType>(3);
		Load_Tracks("exitTracks", types);
		types = static_cast<TrackType>(4);
		Load_Tracks("InactiveLongTracks", types);
		cout << "\t\tNo. of active Short tracks:	= " << activeShortTracks.size() << endl;
		cout << "\t\tNo. of active Long tracks:	= " << activeLongTracks.size() << endl;
		cout << "\t\tNo. of exited tracks:		" << exitTracks.size() << endl;
		cout << "\t\tNo. of inactive tracks:		= " << inactiveTracks.size() << endl;
		cout << "\t\tNo. of inactive Long tracks:	= " << inactiveLongTracks.size() << endl;*/

		if (!(debug_mode == SKIP_PREVIOUS_BACK_STB  && debug_frame_number >= 4 && debug_frame_number <= last)) {
			ConvergencePhase();
		}

		double duration = clock() - start0;

		cout << endl << endl << "\tTotal time taken by STB convergence phase: " << duration / (double) CLOCKS_PER_SEC << "s" << endl;


	}
}

// Initial Phase: for the first 'n = 4' frames, using predictive field for tracking
void STB::InitialPhase(string pfieldfile) {
	int endFrame = last+1;
	if (!(_ipr.triangulationOnly) && !(_ipr.IPROnly)) 
		endFrame = (last >= first + 4) ? first + 4 : last+1;


		for (int frame = first; frame < endFrame; ++frame) {											// ipr on first 4 frames
			std::cout << "IPR on frame " << frame << " of " << last << endl;

//			if (iprFlag || _ipr.triangulationOnly || _ipr.IPROnly) {									// getting 3D positions using ipr (or)
			if (iprFlag && debug_mode != SKIP_IPR) {
				Frame all(_ipr.FindPos3D(imgSequence,frame));
				iprMatched.push_back(all);
//				string IPR_save_path = tiffaddress + "/Tracks/IPRcandidates/" + to_string(frame) +".txt";
//				_ipr.SaveParticlePositions(all.Get_PosDeque(), IPR_save_path);
			}
																
			else {																						// getting 3D positions from .mat files
				string pos3D = "pos3Dframe" + to_string(frame);
				iprMatched.push_back(Load_3Dpoints(pos3D));
			}
		}
																
	if (!(_ipr.triangulationOnly) && !(_ipr.IPROnly)) {												// if the user chooses to run only triangulation or IPR, stop here
																									// else,
		for (int frame = first; frame < endFrame - 1; frame++) {									// using predictive field to get the initial tracks
			int currFrame = frame, nextFrame = frame + 1;
			cout << "STB initial phase tracking b/w frames: " << currFrame << " & " << nextFrame << endl;

			double start = clock();
			PredictiveField pField;
			// getting the predictive velocity field b/w the current and next frames
//			string predict_file_path = tiffaddress + "PredictiveField" + to_string(currFrame - first) + "&"
//												+ to_string(nextFrame - first) + ".txt";
//			if (debug_mode == SKIP_PREDITIVE_FIELD) {
//				pField.Load_field(predict_file_path);
//				if (error == NO_FILE) {
//					cout<<"Predictive file for frame"<<currFrame - first<<"&"<<nextFrame - first<<"can't be opened!"<<endl;
//					goto PredictiveField;
//				}
//			} else {
//				PredictiveField:
				pField.GetPredictiveField(iprMatched[currFrame - first], iprMatched[nextFrame - first], pfieldfile, currFrame);	// either calulating it or from .mat files
//				pField.SaveField(predict_file_path);
//			}
			double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
			cout<<"time to get predictive field:"<< duration << endl;
																									// extending all the tracks that are active in current frame
			start = clock();
			for (deque<Track>::iterator tr = activeShortTracks.begin(); tr != activeShortTracks.end(); ) {
				bool active;
				Position currVelocity(pField.ParticleInterpolation(tr->Last()));
				MakeLink(nextFrame, currVelocity, searchRadiusSTB, *tr, active);

				if (!active) 																		// if the track couldn't find a link, then move the track from active to inactive tracks
					tr = activeShortTracks.erase(tr);
				
				else
					++tr;
			}

			StartTrack(currFrame, pField);															// starting a track for all particles left untracked in current frame
			duration = (clock() - start) / (double) CLOCKS_PER_SEC;
			cout<<"time to link particles in two frames:"<< duration<<endl;
		}
																									// moving all tracks longer than 3 from activeShortTracks to activeLongTracks
		for (deque<Track>::iterator tr = activeShortTracks.begin(); tr != activeShortTracks.end(); ) {
			if (tr->Length() > 3) {
				activeLongTracks.push_back(*tr);
				tr = activeShortTracks.erase(tr);
			}
			else
				++tr;
		}
																									// adding the untracked particles from frame 4 to activeShortTracks
		for (Frame::const_iterator pID = iprMatched[(endFrame - 1) - first].begin(); pID != iprMatched[(endFrame - 1) - first].end(); ++pID) {
			if (!(pID->IsTracked())) {
				Track temp(*pID, endFrame - 1);
				activeShortTracks.push_back(temp);
			}
		}

		cout << "Done with initial phase" << endl << endl;

		cout << "\t\tNo. of active Short tracks:	" << activeShortTracks.size() << endl;
		cout << "\t\tNo. of active Long tracks:	" << activeLongTracks.size() << endl;
		cout << "\t\tNo. of exited tracks:		" << exitTracks.size() << endl;
//		cout << "\t\tNo. of inactive tracks:		" << inactiveTracks.size() << endl;
	}
}

// Convergence phase: after getting tracks of length 4 from initial phase
void STB::ConvergencePhase() {
	int endFrame = last;

	Calibration calib(_ipr.calibfile, ALL_CAMS, _ipr.mindist_2D, _ipr.mindist_3D, ncams);

	for (int currFrame = first + 3; currFrame < endFrame; currFrame++) {							// tracking frame by frame using Wiener / polynomial predictor along with shaking
		string address = tiffaddress + "Tracks/ConvergedTracks/";
		if (debug_mode == SKIP_PREVIOUS_TRACKS && currFrame < debug_frame_number ) {
			if (currFrame < debug_frame_number - 1) continue; // skip those previous frame
				LoadAllTracks(address, to_string(debug_frame_number), 0);

				if (error == NO_FILE) {
					cout<<"The file of tracks can't be opened! The code will load the previous one!\n";
					debug_frame_number --;
					currFrame --;
					error = NONE;
				} else {
					currFrame = debug_frame_number - 1;
					debug_mode = NO_SKIP; //when execute the debug requirement, then label it as done.
				}

		} else {
			// initializing some variables
			int c1 = activeShortTracks.size(), c2 = activeLongTracks.size(), c3 = inactiveTracks.size(), c4 = inactiveLongTracks.size();
			a_as = 0; a_al = 0; a_is = 0; s_as1 = 0; s_as2 = 0; s_as3 = 0; s_as4 = 0; s_al = 0; a_il = 0;

//			for (int i = 0; i < 100; ++i) {
//					linear_fit_error[i] = 0;
//				}

			int nextFrame = currFrame + 1;
			cout << "\tSTB convergence phase tracking at frame: " << nextFrame << endl;
//			cout << "\t\tCheck Points: " << endl;

	// LOADING IMAGES
			deque<string> imgNames;
			for (int i = 0; i < ncams; i++)
				imgNames.push_back(imgSequence[i][nextFrame - 1]);

			Tiff2DFinder t(ncams, _ipr.Get_threshold(), imgNames);
			t.FillPixels(pixels_orig);																	// filled pixels_orig with the px of actual camera images, which will be used for shaking


	// PREDICTION
			cout << "\t\t\tPrediction: ";
			Frame estPos;																				// a frame that saves all the positions estimated from track prediction
			deque<double> estInt;																		// saving the intensity of estimated positions
			auto start = std::chrono::system_clock::now();

			Prediction(nextFrame, estPos, estInt);														// getting the prediction list for all the activLongTrakcs (>4)
			cout<<"Prediction  time:"<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()- start) .count() << "s" << endl;

//			start2 = clock();
//			cout << "Done (" << (clock() - start1)/(double) CLOCKS_PER_SEC << "s)" << endl << "\t\t\tShaking the predictions: " ;
			cout<<"Done.\n";
			/* TESTING */ tempPredictions = estPos;


	// SHAKING THE PREDICTIONS
			start =std::chrono::system_clock::now();
			Shake(estPos, estInt);																		// correcting the predicted positions by shaking, removing wrong / ambiguous predictions and updating the residual images
			cout<<"Shaking  time:"<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start) .count()<< "s" << endl;

	// TESTING PREDICTIONS
			/*_ipr.MatfilePositions(estPos.Get_PosDeque(), tiffaddress + "shakedframe" + to_string(nextFrame));
			_ipr.MatfilePositions(temp.Get_PosDeque(), tiffaddress + "predictionsframe" + to_string(nextFrame));*/
			//if (nextFrame % 100 == 97 || nextFrame % 100 == 98 || nextFrame % 100 == 99 || nextFrame % 100 == 0)
			//	_ipr.MatfileImage(pixels_orig, tiffaddress + "resImg" + to_string(nextFrame));
	// END TESTING

			for (int i = 0; i < estPos.NumParticles(); i++) 											// adding the corrected particle position in nextFrame to its respective track
				activeLongTracks[i].AddNext(estPos[i], nextFrame);

//			start3 = clock();
//			cout << "Done (" << (clock() - start2)/(double) CLOCKS_PER_SEC << "s)" << endl << "\t\t\tShort tracks from residuals: ";
			cout<<"Done.\n";


//			NumDataIO<int> data_io;
//			for (int n = 0; n < ncams; n++) {
//				data_io.SetFilePath("/home/tanshiyong/Documents/Data/SinglePhase/SD78500/Residual_image_code/cam" + to_string(n) + "/frame" + to_string(nextFrame) + ".txt");
//				int* pixel_vector = new int[Npixh * Npixw];
//				for (int i = 0; i < Npixh; i++)
//					for (int j = 0; j < Npixw; j++) {
//						pixel_vector[i * Npixw + j] = pixels_orig[n][j][i];
//					}
//				data_io.SetTotalNumber(Npixw * Npixh);
//				data_io.WriteData(pixel_vector);
//				delete[] pixel_vector;
//			}

//			NumDataIO<int> data_io;
//			for (int n = 0; n < ncams; n++) {
//				data_io.SetFilePath("/home/tanshiyong/Documents/Data/SinglePhase/SD78500//Residual_image_code/cam" + to_string(n+1) + "/frame" + to_string(nextFrame) + ".txt");
//				int* pixel_vector = new int[Npixh * Npixw];
//				for (int i = 0; i < Npixh; i++)
//					for (int j = 0; j < Npixw; j++) {
//						pixel_vector[i * Npixw + j] = pixels_orig[n][i][j];
//					}
//				data_io.SetTotalNumber(Npixw * Npixh);
//				data_io.WriteData(pixel_vector);
//				delete[] pixel_vector;
//			}


	// IPR ON RESIDUALS
			start =std::chrono::system_clock::now();
			_ipr.SetFrameNumber(currFrame + 1);  // frame number is used for debug.
			_ipr.DoingSTB(true);  // to indicate the IPR is called by STB.
			Frame candidates = IPRonResidual(calib, t, pixels_orig, pixels_reproj, pixels_res, estPos);	// applying ipr on residual images to obtain particle candidates
			cout<<"IPR time:"<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start) .count() << "s" << endl;
			cout<<"Total particles:" << candidates.NumParticles()<<endl;

//			string IPR_save_path = tiffaddress + "/Tracks/IPRcandidates/" + to_string(currFrame + 1) +".txt";
//			_ipr.SaveParticlePositions(candidates.Get_PosDeque(), IPR_save_path);

//			NumDataIO<double> data_io;
//			double* pos_vec = new double[candidates.NumParticles() * 3];
//			for (int i = 0; i < candidates.NumParticles(); ++i) {
//				pos_vec[i * 3] = candidates[i].X();
//				pos_vec[i * 3 + 1] = candidates[i].Y();
//				pos_vec[i * 3 + 2] = candidates[i].Z();
//			}
//			data_io.SetFilePath("/home/tanshiyong/Documents/Data/Single-Phase/SD50000/frame" + to_string(nextFrame) + ".txt");
//			data_io.SetTotalNumber(candidates.NumParticles() * 3);
//			data_io.WriteData(pos_vec);

	// TESTING IPR AND TRACKING FROM RESIDUALS
//			if (nextFrame % 100 == 97 || nextFrame % 100 == 98 || nextFrame % 100 == 99 || nextFrame % 100 == 0)
//				_ipr.SaveParticlePositions(candidates.Get_PosDeque(), tiffaddress + "pos3Dresframe" + to_string(nextFrame));
	// END TESTING
			start  = std::chrono::system_clock::now();

			/*
			 *  Modified by Shiyong Tan
			 *  Parallelizing this section in order to increase the processing speed
			 *  Introduce variable erase_vector, candidate_used to label the tracks and candidate points which should be deleted after linking.
			 *  It is used to avoid to erase them in the loop which would result in change of size of the loop that hinders parallelization.
			 *  Start:
			 */
																										// trying to link each activeShortTrack with a particle candidate
			unsigned int num_shorttrack = activeShortTracks.size();
			bool* erase_vector = new bool[num_shorttrack];  				// Vector to label the tracks to be deleted.
			unsigned int num_candidate = candidates.NumParticles();
			bool* candidate_used = new bool[num_candidate];			// Vector to label the candidate points that are used and need to be deleted.
			for (unsigned int i = 0; i < num_candidate; ++i) candidate_used[i] = false;
#pragma omp parallel
			{
#pragma omp for
			for (int i = 0; i < num_shorttrack; ++i) {
				deque<Track>::iterator tr = activeShortTracks.begin() + i;
				bool erase = false;
				MakeShortLinkResidual(nextFrame, candidates, tr, 5, &erase, candidate_used);  // erase and candidate_used are updated for every loop
				erase_vector[i] = erase;
			}
			}

			// Delete labeled tracks
			int shift = 0;
			for (int i = 0; i < num_shorttrack; ++i) {
				deque<Track>::iterator tr = activeShortTracks.begin();
				if (erase_vector[i]) {
					activeShortTracks.erase(tr + i - shift);
					++ shift;
				}
			}
			delete[] erase_vector;
			// Delete labeled candidates
			shift = 0;
			for (int i = 0; i < num_candidate; ++i) {
				if (candidate_used[i]) {
					candidates.Delete(i - shift);
					++ shift;
				}
			}
			delete[] candidate_used;
//			for (deque<Track>::iterator tr = activeShortTracks.begin(); tr != activeShortTracks.end(); )
//				MakeShortLinkResidual(nextFrame, candidates, tr, 5);
			// END

			cout<<"linking short track time:"<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()- start) .count() << "s" << endl;
			start = std::chrono::system_clock::now();
																										// moving all activeShortTracks longer than 3 particles to activeLongTracks
			for (deque<Track>::iterator tr = activeShortTracks.begin(); tr != activeShortTracks.end(); ) {

				if (tr->Length() > 3) {
					activeLongTracks.push_back(*tr);
					tr = activeShortTracks.erase(tr);
					s_as4++; a_al++;
				}
				else
					++tr;
			}
//			start4 = clock();
//			cout << "Done (" << (clock() - start3) / (double) CLOCKS_PER_SEC << "s)" << endl << "\t\t\tPruning the tracks: ";
			cout<<"Done.\n";

	// PRUNING / ARRANGING THE TRACKS
			double thresh = 1.5 * largestShift;
//			int p = 0;
			int r1 =0, r2 = 0, r3 = 0, r4 = 0;
			unsigned int num_longtrack = activeLongTracks.size();
			erase_vector  = new bool[num_longtrack];
			for (unsigned int i = 0; i < num_longtrack; ++i) erase_vector[i] = false;

//			start = std::chrono::system_clock::now();
////			unique_ptr<double> dist_map (new double(num_longtrack * num_longtrack));
//			double* dist_map = new double[num_longtrack * num_longtrack];
//			double search_radius = 5;
//			for (unsigned int i = 0; i < num_longtrack; ++ i) {
//				Position last = activeLongTracks[i].Last();
//				double x_l = last.X() - search_radius, x_u = last.X() + search_radius;
//				double y_l = last.Y() - search_radius, y_u = last.Y() + search_radius;
//				double z_l = last.Z() - search_radius, z_u = last.Z() + search_radius;
//
//			#pragma omp parallel //num_threads(8)
//									{
//			#pragma omp for
//				for (unsigned int j = i; j < num_longtrack; ++j) {
//					dist_map[i * num_longtrack + j] = -1;
//					dist_map[j * num_longtrack + i] = -1;
//					double x_test = activeLongTracks[j].Last().X();
//
//			//		double d = pow(Distance(tr->Last(), activeLongTracks[i].Last()), 0.5);
//			//		if (fabs(tr.Last().X() - activeLongTracks[i].Last().X()) < search_radius)
//			//			if (fabs(tr.Last().Y() - activeLongTracks[i].Last().Y()) < search_radius)
//			//				if (fabs(tr.Last().Z() - activeLongTracks[i].Last().Z()) < search_radius) {
//					if (x_test > x_l && x_test < x_u) {
//						double y_test = activeLongTracks[j].Last().Y();
//							if (y_test > y_l && y_test < y_u) {
//								double z_test = activeLongTracks[j].Last().Z();
//								if (z_test > z_l && z_test < z_u ) {
//									double d = pow(Distance(last, activeLongTracks[i].Last()), 0.5);
//									dist_map[i * num_longtrack + j] = d;
//									dist_map[j * num_longtrack + i] = d;
//								}
//							}
//					}
//				}
//									}
//			}
//			cout<<"Calculating dist map:"<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()- start) .count() << "s" << endl;

#pragma omp parallel
			{
#pragma omp for
			for (unsigned int i = 0; i < num_longtrack; ++ i) {
//				double d1 = pow(Distance(activeLongTracks[i].Last(), activeLongTracks[i].Penultimate()),0.5),
//						d2 = pow(Distance(activeLongTracks[i].Penultimate(), activeLongTracks[i].Antepenultimate()),0.5);
//				double threshRel = maxRelShiftChange * d2; //, length = activeLongTracks[i].Length();
//				if (d1 > thresh) {
//					erase_vector[i] = true;
//					++ r1;
//				}
//				else if (fabs(d1 - d2) > maxAbsShiftChange || fabs(d1 - d2) > threshRel) {
//					erase_vector[i] = true;
//					++ r2;
//				}
//				if (!CheckVelocity(activeLongTracks[i])) {
//					erase_vector[i] = true;
//					++ r1;
//				}
				 if (!CheckLinearFit(activeLongTracks[i])) {
					++ r3;
					erase_vector[i] = true;
				}
			}
			}

			shift = 0;
			for (unsigned int i = 0; i < num_longtrack; ++ i) {
				deque<Track>::iterator tr = activeLongTracks.begin();
				if (erase_vector[i]) {
					if (activeLongTracks[i - shift].Length() >= 7) {
						inactiveLongTracks.push_back(activeLongTracks[i]);
						activeLongTracks.erase(tr + i - shift);
						s_al++; a_il++;
						++ shift;
					} else {
						activeLongTracks.erase(tr + i - shift);
						s_al++; a_is++;
						++ shift;
					}
				}
			}
			delete[] erase_vector;

//			for (deque<Track>::iterator tr = activeLongTracks.begin(); tr != activeLongTracks.end(); ) {
////				p ++;
////				cout<<p<<endl;
//				double d1 = pow(Distance(tr->Last(), tr->Penultimate()),0.5), d2 = pow(Distance(tr->Penultimate(), tr->Antepenultimate()),0.5);
//				double threshRel = maxRelShiftChange*d2, length = tr->Length();
//																										// moving all activeLongTracks with displacement more than LargestExp shift to inactiveTracks
//				if (d1 > thresh) {
////					inactiveTracks.push_back(*tr);
//					tr->DeleteBack();
//					if (length >= 7) {
//						inactiveLongTracks.push_back(*tr);
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_il++;
//					} else {
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_is++;
//					}
//					r1 ++;
//					continue;
//				}
//																										// moving all activeLongTracks with large change in particle shift to inactive/inactiveLong tracks
//				else if (fabs(d1 - d2) > maxAbsShiftChange || fabs(d1 - d2) > threshRel) {
//					tr->DeleteBack();
//					if (length >= 7) {
//						inactiveLongTracks.push_back(*tr);
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_il++;
//					}
//					else {
////						inactiveTracks.push_back(*tr);
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_is++;
//					}
//					r2 ++;
//					continue;
//				}
//
//				if (!CheckVelocity(tr)) { // check velocity
//					if (length >= 7) {
//						inactiveLongTracks.push_back(*tr);
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_il++;
//					} else {
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_is++;
//					}
//					r4 ++;
//					continue;
//				}
////				else if (!CheckAcceleration(tr)) {
////						tr->DeleteBack();
////						if (length >= 7) {
////							inactiveLongTracks.push_back(*tr);
////						}
////						tr = activeLongTracks.erase(tr);
////						s_al++; a_is++;
////					}
//
//				else if (!CheckLinearFit(tr)) {
//					if (length >= 7) {
//						inactiveLongTracks.push_back(*tr);
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_il++;
//					} else {
//						tr = activeLongTracks.erase(tr);
//						s_al++; a_is++;
//					}
//					r3 ++;
//					continue;
//				}
//
//				else
//					++tr;
//			}
//			GetAccThred(); // Get the threshold for acceleration check


			for (Frame::const_iterator pID = candidates.begin(); pID != candidates.end(); ++pID) {		// adding all the untracked candidates to a new short track
				Track startTrack(*pID, nextFrame);
				activeShortTracks.push_back(startTrack); a_as++;
			}
			cout<<"Pruning time:"<< std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start).count() << "s" << endl;

//			cout << "Done (" << (clock() - start4) / (double) CLOCKS_PER_SEC << "s)" << endl;

			cout << "\t\tNo. of active Short tracks:	" << c1 << " + " << a_as << " - (" << s_as1 << " + " << s_as2 << " + " << s_as3 << " + " << s_as4 << ") = " << activeShortTracks.size() << endl;
			cout << "\t\tNo. of active Long tracks:	" << c2 << " + " << a_al << " - " << s_al << " = " << activeLongTracks.size() << endl;
			cout << r1 << "+" << r2 << "+" << r3 << endl;
			cout << "\t\tNo. of exited tracks:		 = " << exitTracks.size() << endl;
//			cout << "\t\tNo. of inactive tracks:		" << c3 << " + " << a_is << " = " << inactiveTracks.size() << endl;
			cout << "\t\tNo. of inactive Long tracks:	" << c4 << " + " << a_il << " = " << inactiveLongTracks.size() << endl;

//			cout << "\t\tTime taken for STB at frame " << nextFrame << ": " << (clock() - start0) / (double) CLOCKS_PER_SEC << "s" << endl;



			if (to_save_data) {
				MatTracksSave(address, to_string(nextFrame), 0);
			} else {
				//time_t t = time(0);
				if (nextFrame % 100 == 0 || nextFrame == endFrame) {  // to debug, every frame should be saved
					cout << "\tSaving the tracks" << endl;

				MatTracksSave(address, to_string(nextFrame), 0);
				if (nextFrame % 100 == 0) {
					// remove previous files except anyone that is multiple of 500 for active long track and exit track
					std::remove((address + "ActiveLongTracks" + to_string(nextFrame - 100) + ".txt").c_str());
					std::remove( (address + "ActiveShortTracks" + to_string(nextFrame - 100) + ".txt").c_str());
					// remove previous files except anyone that is multiple of 500 for active long track and exit track
					if ((nextFrame - 100) % 500 != 0) {
						std::remove( (address + "InactiveLongTracks" + to_string(nextFrame - 100) + ".txt").c_str());
						std::remove( (address + "ExitTracks" + to_string(nextFrame - 100) + ".txt").c_str());
					}
				} else if (nextFrame == endFrame && nextFrame % 100 != 0) {
					std::remove((address + "ActiveLongTracks" + to_string((int)(nextFrame / 100) * 100) + ".txt").c_str());
					std::remove( (address + "ActiveShortTracks" + to_string((int)(nextFrame / 100) * 100) + ".txt").c_str());
					if (((int)(nextFrame / 100) * 100) % 500 != 0) {
						std::remove( (address + "InactiveLongTracks" + to_string((int)(nextFrame / 100) * 100) + ".txt").c_str());
						std::remove( (address + "ExitTracks" + to_string((int)(nextFrame / 100) * 100) + ".txt").c_str());
					}
				}


				}
			}

			/*
			 * Modified by Shiyong Tan
			 * Save the inacitve long tracks for every 500 frames and empty it to avoid endless expansion of this variable
			 * Start:
			 */
			if (nextFrame % 500 == 0) {
				// save the inactive long tracks
				string s = to_string(nextFrame);
				string X5 = "InactiveLongTracks" + s;
				SaveTrackToTXT(inactiveLongTracks, address + X5);
				// empty inactiveLongTracks
				if (nextFrame != endFrame) inactiveLongTracks.erase(inactiveLongTracks.begin(), inactiveLongTracks.end());
				// save the inactive long tracks
				string X6 = "ExitTracks" + s;
				SaveTrackToTXT(exitTracks, address + X6);
				// empty inactiveLongTracks
				if (nextFrame != endFrame) exitTracks.erase(exitTracks.begin(), exitTracks.end());
			}
		}
	}
}

// a function to start a track
void STB::StartTrack(int frame, PredictiveField& pField) {
	for (Frame::const_iterator pID = iprMatched[frame - first].begin(); pID != iprMatched[frame - first].end(); ++pID) {																									
		if (!(pID->IsTracked())) {																// if the particle is not a part of a track
			Track initialTrack;																		// then start a track
			
			iprMatched[frame - first][pID.where()].SetTracked();									// setting the particle as tracked and adding it to a track
			initialTrack.AddNext(*pID, frame);
			//cout << "(X,Y,Z) = (" << pID->X() << "," << pID->Y() << "," << pID->Z() << ")" << " and track is active: " << initialTrack.IsActive() << endl;
			
			Position currVelocity(pField.ParticleInterpolation(*pID));								// interpolating the velocity field to pID
			
			bool activeTrack;
			
			MakeLink(frame + 1, currVelocity, searchRadiusSTB, initialTrack, activeTrack);			// finding a link for pID in next frame
			
			if (activeTrack)
				activeShortTracks.push_back(initialTrack);

			// setting the track as active: if it finds an element in next frame / inactive: if it doesn't 
			//allTracks.back().SetActive(activeTrack);
		}	
	}
}

// a function to make a link using predictive velocity field
void STB::MakeLink(int nextFrame, const Position& currVelocity, double radius, Track& track, bool& activeTrack)
{
		
		Position estimate;																			// we'll need an estimated future position
		
		estimate = track.Last() + currVelocity;														// the estimated point	
		
		int len = track.Length();																	// length of the track
		
		pair<Frame::const_iterator, float> cost;
		
		if (len == 1) {																		// if the track has only one particle
																									// choosing the particle closest to the estimate
			cost = ComputeCost(iprMatched[nextFrame - first], iprMatched[nextFrame - first], radius, estimate, currVelocity, currVelocity, true);
		}
		else if (len >= 2) {																// if it has more than 2 particles		
			Position prevVelocity = track.Last() - track.Penultimate();								// using the velocity b/w previous 2 frames
			/*if (nextFrame - first < iprMatched.size() - 1) 
				cost = ComputeCost(iprMatched[nextFrame - first], iprMatched[nextFrame - first + 1], radius, estimate, prevVelocity, currVelocity, false);
			
			else */
				cost = ComputeCost(iprMatched[nextFrame - first], iprMatched[nextFrame - first], radius, estimate, prevVelocity, currVelocity, true);
		}
						
		if (cost.second == UNLINKED) {														// if no link found in next frame		
			//track.SetActive(false);																// set as inactive track
			activeTrack = false;
		}
		else {																				// if a link is found																									
			iprMatched[nextFrame - first][cost.first.where()].SetTracked();							// setting the particle from nextFrame as tracked		
			track.AddNext(*cost.first, nextFrame);													// adding the particle to current track
			activeTrack = true;																		// setting the track as still active
			//cout << "(X,Y,Z) = (" << cost.first->X() << "," << cost.first->Y() << "," << cost.first->Z() << ")" << " and track is active: " << track.IsActive() << endl;
		}
}

pair<Frame::const_iterator, float> STB::ComputeCost(Frame& fr1, Frame& fr2, double radius,
									const Position& estimate,
									const Position& prevVelocity,
									const Position& currVelocity,
									bool stopflag,
									bool* candidate_used)
{
																									// find possible continuations in the next frame.

	deque<Frame::const_iterator> matches;															// find possible matches
	deque<float> match_costs;
	deque<unsigned int>match_index;
	float mincost = radius * radius;

	unsigned int fr_no = 0;
	Frame::const_iterator fitend = fr1.end();
	for (Frame::const_iterator fit = fr1.begin(); fit != fitend; ++fit, ++fr_no) {
		if (candidate_used[fr_no])  {
//			++fr_no;
			continue;
		}
//		++fr_no;
		float mag = Distance(estimate, *fit);
		if (mag > mincost) {
			continue;
		}
		matches.push_back(fit);
		match_index.push_back(fr_no);   //Shiyong Tan:  get the index of the candidate in order to label them later.
		if (stopflag == true) {
			match_costs.push_back(mag);																// don't look into the future any more
			if (mincost > mag) {
				mincost = mag;
			}
		}
		else {
			Position acceleration = currVelocity - prevVelocity;									// project again!
			Position new_estimate = *fit + currVelocity + 0.5 * acceleration;
			pair<Frame::const_iterator, float> cost = ComputeCost(fr2, fr2, radius, new_estimate, prevVelocity, currVelocity, true, candidate_used);
			match_costs.push_back(cost.second);
			if (mincost > cost.second) {
				mincost = cost.second;
			}
		}
	}

	deque<Frame::const_iterator>::const_iterator m_end = matches.end();
	pair<Frame::const_iterator, float> retval = make_pair(fitend-1, -1);
	deque<float>::const_iterator mc = match_costs.begin();
	unsigned int index_used = 0;
	for (deque<Frame::const_iterator>::const_iterator m = matches.begin(); m != m_end; ++m, ++mc, ++index_used) {
		if (*mc > mincost) {
			continue;
		}
		retval = make_pair(*m, *mc);
		candidate_used[match_index[index_used]]= true;
	}

	return retval;
}

pair<Frame::const_iterator, float> STB::ComputeCost(Frame& fr1, Frame& fr2, double radius,
									const Position& estimate,
									const Position& prevVelocity,
									const Position& currVelocity,
									bool stopflag)
{
																									// find possible continuations in the next frame.

	deque<Frame::const_iterator> matches;															// find possible matches
	deque<float> match_costs;
	float mincost = radius * radius;

	Frame::const_iterator fitend = fr1.end();
	for (Frame::const_iterator fit = fr1.begin(); fit != fitend; ++fit) {
		float mag = Distance(estimate, *fit);
		if (mag > mincost) {
			continue;
		}
		matches.push_back(fit);
		if (stopflag == true) {
			match_costs.push_back(mag);																// don't look into the future any more
			if (mincost > mag) {
				mincost = mag;
			}
		}
		else {			
			Position acceleration = currVelocity - prevVelocity;									// project again!
			Position new_estimate = *fit + currVelocity + 0.5 * acceleration;
			pair<Frame::const_iterator, float> cost = ComputeCost(fr2, fr2, radius, new_estimate, prevVelocity, currVelocity, true);
			match_costs.push_back(cost.second);
			if (mincost > cost.second) {
				mincost = cost.second;
			}
		}
	}

	deque<Frame::const_iterator>::const_iterator m_end = matches.end();
	pair<Frame::const_iterator, float> retval = make_pair(fitend-1, -1);
	deque<float>::const_iterator mc = match_costs.begin();
	for (deque<Frame::const_iterator>::const_iterator m = matches.begin(); m != m_end; ++m, ++mc) {
		if (*mc > mincost) {
			continue;
		}
		retval = make_pair(*m, *mc);
	}

	return retval;
}

void STB::Prediction(int frame, Frame& estPos, deque<double>& estInt) {
	int t = frame;																			// time
	vector<vector<double>> predCoeff(3);													// using polynomial fit / Wiener filter to predict the particle position at nextFrame
	vector<string> direction = { "X", "Y", "Z" };
	vector<double> est(3);
//	deque<Track>::const_iterator tr_end = activeLongTracks.end();								// for each track in activeLongTracks (at least 4 particles long)
	int debug_index = 0;
	int debug_index1 = 0;
	/*
	 * Modified by: Shiyong Tan
	 * After boundary Check, if the prediction is out of boundary, then move the track into exit track.
	 * start:
	 */
	for (deque<Track>::iterator tr = activeLongTracks.begin(); tr != activeLongTracks.end(); ) {
		++ debug_index;
		for (int i = 0; i < 3; i++) {															// getting predictor coefficients for X->0, Y->1 and Z->2
			if (tr->Length() < 4) {																// if length is 4 or 5, use all points to get 2nd degree polynomial
//				predCoeff[i] = Polyfit(*tr, direction[i], tr->Length(), 2);
//				est[i] = predCoeff[i][0] + predCoeff[i][1] * t + predCoeff[i][2] * pow(t, 2);
				est[i] = LMSWienerPredictor(*tr, direction[i], 3);
			}
			else if (tr->Length() < 6) {														// if length is 6 to 10, use all points to get 3rd degree polynomial
//				predCoeff[i] = Polyfit(*tr, direction[i], tr->Length(), 3);
//				//cout << "break" << tr->Length() << endl;
//				est[i] = predCoeff[i][0] + predCoeff[i][1] * t + predCoeff[i][2] * pow(t, 2) + predCoeff[i][3] * pow(t, 3);
				est[i] = LMSWienerPredictor(*tr, direction[i], tr->Length() - 1);
			}
			else {																				// if length is more than 11, use last 10 points to get 3rd degree polynomial
//				predCoeff[i] = Polyfit(*tr, direction[i], 10, 3);
//				est[i] = predCoeff[i][0] + predCoeff[i][1] * t + predCoeff[i][2] * pow(t, 2) + predCoeff[i][3] * pow(t, 3);
				est[i] = LMSWienerPredictor(*tr, direction[i], 5);
			}
		}
		Position estimate(est[0], est[1], est[2]);												// estimated position at nextFrame
		if (boundary_check.Check(estimate)) { // if the estimate particle is inside the boundary
			estInt.push_back(1);																	// setting initial intensity to 1
			estPos.Add(estimate);
			++tr;
			++ debug_index1;
		} else { //the track should be put into exit tracks since it is outside the boundary
			exitTracks.push_back(*tr);
			tr = activeLongTracks.erase(tr);
			s_al ++;
		}
		//End
	}
}

// convergence phase
// filter for prediction 
// polynomial fit
vector<double> STB::Polyfit(Track tracks, string direction, int datapoints, int polydegree) {

	int size = tracks.Length();

	double* t = new double[datapoints];
	double* x = new double[datapoints];

	for (int i = 0; i < datapoints; i++) {
		t[i] = tracks.GetTime(size - i - 1);	//time

		if (direction == "X" || direction == "x")
			x[i] = tracks[size - i - 1].X();	//x-values or
		else if (direction == "Y" || direction == "y")
			x[i] = tracks[size - i - 1].Y();	//y-values or
		else if (direction == "Z" || direction == "z")
			x[i] = tracks[size - i - 1].Z();	//z-values
	}

	// n = 3 is the degree;
	double* T = new double[2 * polydegree + 1];                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	for (int i = 0; i < 2 * polydegree + 1; i++) {
		T[i] = 0;
		for (int j = 0; j<datapoints; j++)
			T[i] = T[i] + pow(t[j], i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	}
									            
	double** B = new double*[polydegree + 1];	//B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
	for (int i = 0; i < polydegree + 1; i++)
		B[i] = new double[polydegree + 2];
	double* a = new double[polydegree + 1];
	for (int i = 0; i < polydegree + 1; i++)
		a[i] = 0;

	for (int i = 0; i < polydegree + 1; i++) 
		for (int j = 0; j < polydegree + 1; j++)
			B[i][j] = T[i + j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix

	double* X = new double[polydegree + 1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	for (int i = 0; i < polydegree + 1; i++) {
		X[i] = 0;
		for (int j = 0; j<datapoints; j++)
			X[i] = X[i] + pow(t[j], i)*x[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	}

	for (int i = 0; i < polydegree + 1; i++)
		B[i][polydegree + 1] = X[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)

	polydegree = polydegree + 1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations
	//cout << "\nThe Normal(Augmented Matrix) is as follows:\n";

	for (int i = 0; i < polydegree; i++)                     //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
		for (int k = i + 1; k < polydegree; k++) 
			if (B[i][i] < B[k][i]) 
				for (int j = 0; j <= polydegree; j++) {
					double temp = B[i][j];
					B[i][j] = B[k][j];
					B[k][j] = temp;
				}			

	for (int i = 0; i < polydegree - 1; i++)             //loop to perform the gauss elimination
		for (int k = i + 1; k < polydegree; k++)	{
			double s = B[k][i] / B[i][i];
			for (int j = 0; j <= polydegree; j++)
				B[k][j] = B[k][j] - s*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables		
		}
	
	for (int i = polydegree - 1; i >= 0; i--)                //back-substitution
	{                        //x is an array whose values correspond to the values of x,y,z..
		a[i] = B[i][polydegree];                //make the variable to be calculated equal to the rhs of the last equation
		for (int j = 0; j < polydegree; j++) {
			if (j != i)            //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
				a[i] = a[i] - B[i][j] * a[j];
		}
		a[i] = a[i] / B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
	}

	vector<double> coeff;
	for (int i = 0; i < polydegree; i++) 
		coeff.push_back(a[i]);
	
	polydegree = polydegree - 1;
	for (int i = 0; i < polydegree + 1; i++)
		delete[] B[i];

	delete[] T, x, t, X,a, B;

	return coeff;
}

// predict the next point with Wiener Predictor using LMS algorithm
double STB::LMSWienerPredictor(Track tracks, string direction, int order) {
	int size = tracks.Length();
	double* series = new double[order + 1];

	// get the series to be predicted
	for (unsigned int i = 0; i < order + 1; ++i) {
		if (direction == "X" || direction == "x")
			series[i] = tracks[size - order - 1 + i].X();	//x-values or
		else if (direction == "Y" || direction == "y")
			series[i] = tracks[size - order - 1 + i].Y();	//y-values or
		else if (direction == "Z" || direction == "z")
			series[i] = tracks[size - order - 1 + i].Z();	//z-values
	}

	// Wiener filter does badly near zero
	bool shift_label = false;
	double shift = 10;
	if (fabs(series[order]) < 1) {
		// make a shift to avoid zero-plane prediction
		shift_label = true;
		for (unsigned int i = 0; i < order + 1; ++i) {
			series[i] = series[i] + shift;
		}
	}

	double* filter_param = new double[order]; // the filter parameter
	for (unsigned int i = 0; i < order; ++i) filter_param[i] = 0; // initialize the filter
	// calculate  the step
	double sum = 0;
	for (unsigned int i = 0; i < order; ++i) {
		sum = sum + series[i] * series[i];
	}
	double step = 1 / sum;

	double prediction = 0;
	for (unsigned int i = 0; i < order; ++i) {
		prediction = prediction + filter_param[i] * series[i];
	}

	double error = series[order] - prediction;

	while (fabs(error) > 1e-8) {
		for (unsigned int i = 0; i < order; ++i) {
			filter_param[i] = filter_param[i] + step * series[i] * error;
		}
		//calculate the prediction using the new filter parameters
		prediction = 0;
		for (unsigned int i = 0; i < order; ++i) {
			prediction = prediction + filter_param[i] * series[i];
		}
		error = series[order] - prediction;
	}
	prediction = 0;
	for (unsigned int i = 0; i < order; ++i) {
		prediction = prediction + filter_param[i] * series[i + 1];
	}
	if (shift_label) {
		prediction = prediction - shift;
	}
	delete[] series;
	delete[] filter_param;
	return prediction;
}

// For shaking the predicted estimates
void STB::Shake(Frame& estimate, deque<double>& intensity) {
	
	if (estimate.NumParticles() > 0) {
		for (int index = 0; index < estimate.NumParticles(); index++) {						// adding 2D image centers and intensity data to the estimates 
			_ipr.FullData(estimate[index], intensity[index], ncams, ALL_CAMS);
		}

		deque<int> ignoreCam = Rem(estimate, intensity, _ipr.mindist_3D);					// removing ambiguous, ghost and particles that disappear on at least 2 cams

		for (int loopInner = 0; loopInner < _ipr.it_innerloop; loopInner++) {
			double del;
//			if (loopInner < 1)  del = config.shaking_shift;			// initial shake TODO
//			else if (loopInner < 5)  del = config.shaking_shift / pow(2,loopInner - 1);//_ipr.mindist_2D/10;	// normal shakes TODO
//			else  del = config.shaking_shift/100;

//			double mindist_2D = _ipr.Get_mindist2D();

			if (loopInner < 1)  del = config.shaking_shift;
			else if (loopInner < 5)  del = config.shaking_shift / pow(2,loopInner - 1);
			else  del = config.shaking_shift / 20;

			_ipr.ReprojImage(estimate, OTFcalib, pixels_reproj, STBflag);
//			_ipr.ReprojImage(estimate, OTFcalib, pixels_reproj, 0.5);					// adding the estimated particles to the reprojected image

			for (int n = 0; n < ncams; n++) 												// updating the residual image by removing the estimates
				for (int i = 0; i < Npixh; i++)
					for (int j = 0; j < Npixw; j++) {
						int residual = pixels_orig[n][i][j] - abs(pixels_reproj[n][i][j]);
						pixels_res[n][i][j] =  residual; //(residual < 0) ? 0 : residual;
					}

//			// output the residual image
//			NumDataIO<int> image_output;
//			string save_path;
//			for (int n = 0; n < ncams; n++) {
//				save_path = "/home/sut210/Documents/Experiment/EXP6/resimg" + to_string(n) +".txt";
//				image_output.SetFilePath(save_path);
//				image_output.SetTotalNumber(Npixh * Npixw);
//				int* pixel_array = new int[Npixh * Npixw];
//				for (int i = 0; i < Npixh; i++)
//					for (int j = 0; j < Npixw; j++)
//						pixel_array[i * Npixw + j] = pixels_res[n][i][j];
//				image_output.WriteData(pixel_array);
//			}


//			int index = 0;																	// correcting the estimated positions and their intensity by shaking
//			double start = clock();
#pragma omp parallel num_threads(24)
						{
//							int TID = omp_get_thread_num();
//							cout<<"Thread %d is runing\n"<<TID<<endl;
#pragma omp for
//			for (Frame::const_iterator pID = estimate.begin(); pID != estimate.end(); ++pID) {
			for (int i = 0; i < estimate.NumParticles(); ++i) {
				Frame::const_iterator pID = estimate.begin() + i;
				OTF otf_calib(OTFcalib);
				Shaking s(ncams, ignoreCam[i], otf_calib, Npixw, Npixh, _ipr.psize, del, *pID, cams, pixels_res, intensity[i]);
				estimate[i] = s.Get_posnew();
				intensity[i] = s.Get_int();
//				index++;
			}
						}
//			double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
//			cout<<"time to do shaking for each loop:"<<duration<<endl;
		}

		for (int index = 0; index < estimate.NumParticles(); index++) {						// adding 2D image centers and intensity data after shaking
			_ipr.FullData(estimate[index], intensity[index], ncams, ALL_CAMS);
		}

		Rem(estimate, intensity, _ipr.mindist_3D);											// removing ambiguous particles and particles that did not find a match on the actual images

		_ipr.ReprojImage(estimate, OTFcalib, pixels_reproj, STBflag);						// updating the reprojected image
//		_ipr.ReprojImage(estimate, OTFcalib, pixels_reproj, 1.5);

		for (int n = 0; n < ncams; n++) 													// updating the residual image
			for (int i = 0; i < Npixh; i++)
				for (int j = 0; j < Npixw; j++) {
					int residual = (pixels_orig[n][i][j] - fpt*pixels_reproj[n][i][j]);		// using the multiplying factor (fpt) to remove all traces of tracked particles
					pixels_orig[n][i][j] = (residual < 0) ? 0 : residual;
				}


	}
	
}

// removing ghost,ambiguous and particles leaving measurement domain
deque<int> STB::Rem(Frame& pos3D, deque<double>& int3D, double mindist_3D) {
	double thresh3D = 4 * mindist_3D * mindist_3D;
	double threshShift = pow(largestShift, 2);
	for (int i = 0; i < pos3D.NumParticles(); i++) {									// deleting particles that are very close to each other
		for (int j = i + 1; j < pos3D.NumParticles(); ) {
			if (Distance(pos3D[i], pos3D[j]) < thresh3D) {
				pos3D.Delete(j); int3D.erase(int3D.begin() + j); tempPredictions.Delete(j);
//				inactiveTracks.push_back(activeLongTracks[j]); TODO					// shifting the corresponding activeLongTrack to inactiveTracks
				if (activeLongTracks.at(j).Length() >= 7) {
					inactiveLongTracks.push_back(activeLongTracks.at(j));
					activeLongTracks.erase(activeLongTracks.begin() + j);
					s_al++; a_il++;
				} else {
					activeLongTracks.erase(activeLongTracks.begin() + j);
					s_al++; a_is++;
				}
			}
			else
				j++;
		}
	}

	int ghost = 0;																		
	double avgIntTemp = 0, avgInt = 0, count = 0;										// deleting based on intensity
	for (int i = 0; i < int3D.size(); i++) 
		if (int3D[i] > 0) 
			avgIntTemp = avgIntTemp + int3D[i];
		
	avgInt = avgIntTemp;
	avgIntTemp = avgIntTemp / int3D.size();

	for (int i = 0; i < int3D.size(); i++)												// removing the outliers (very bright particles and very dull particles for mean)
		if (int3D[i] > 30*avgIntTemp) {
			avgInt = avgInt - int3D[i];
			count++;
		}
	avgInt = avgInt / (int3D.size()-count);

	for (int index = int3D.size() - 1; index >= 0; index--) {							// removing a particle if its intensity falls below a % of the mean intensity
		if (int3D[index] < intensityLower * avgInt) {
			//ghost3D.push_back(pos3D[index]); 
			pos3D.Delete(index); int3D.erase(int3D.begin() + index); tempPredictions.Delete(index);
																						// shifting the corresponding activeLongTrack to inactiveLong or inactiveTracks
			if (activeLongTracks[index].Length() >= 7) { // && Distance(activeLongTracks[index].Last(), activeLongTracks[index].Penultimate()) < threshShift) {
				inactiveLongTracks.push_back(activeLongTracks[index]);
				a_il++;
			}
			else {
//				inactiveTracks.push_back(activeLongTracks[index]);						// shifting the corresponding activeLongTrack to inactiveTracks
				a_is++;
			}
			activeLongTracks.erase(activeLongTracks.begin() + index);
			 s_al++;
			ghost++;
		}
	}

	deque<int> ignoreCam(pos3D.NumParticles());
	double intensityThresh = 0.25*_ipr.Get_threshold();
	// new version: Adding a check of intensity 
		// For a given predicted particles, if the intensity on one original image falls below a quarter of peak intensity threshold, ignore that camera for shaking
		// If it falls below on two cameras, move the track to inactive or inactiveLongTracks

	for (int i = 0; i < pos3D.NumParticles(); ) {										// deleting particles that are outside the image bounds on atleast 2 cameras
		double xlim = ((double)(Npixw - 1)) / 2, ylim = ((double)(Npixh - 1)) / 2;		// if the particle disappears only on 1 cam, ignore that cam and perform shaking
		int leftcams = 0, belowIntensity = 0; double shiftThreshold = pow(largestShift, 2);					
		ignoreCam[i] = ALL_CAMS;

		// out of bounds or below intensity threshold
		if (abs(pos3D[i].X1() - xlim) > xlim || abs(pos3D[i].Y1() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 0;
		}
//		else if (pixels_orig[0][(int)round(pos3D[i].Y1())][(int)round(pos3D[i].X1())] < intensityThresh) {
//			belowIntensity++; ignoreCam[i] = 0;
//		}

		if (abs(pos3D[i].X2() - xlim) > xlim || abs(pos3D[i].Y2() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 1;
		}
//		else if (pixels_orig[1][(int)round(pos3D[i].Y2())][(int)round(pos3D[i].X2())] < intensityThresh) {
//			belowIntensity++; ignoreCam[i] = 1;
//		}

		if (abs(pos3D[i].X3() - xlim) > xlim || abs(pos3D[i].Y3() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 2;
		}
//		else if (pixels_orig[2][(int)round(pos3D[i].Y3())][(int)round(pos3D[i].X3())] < intensityThresh) {
//			belowIntensity++; ignoreCam[i] = 2;
//		}

		if (abs(pos3D[i].X4() - xlim) > xlim || abs(pos3D[i].Y4() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 3;
		}
//		else if (pixels_orig[3][(int)round(pos3D[i].Y4())][(int)round(pos3D[i].X4())] < intensityThresh) {
//			belowIntensity++; ignoreCam[i] = 3;
//		}

		if (leftcams >= 2) {											// if the particle disappears on at least 2 cams (out of bounds only)
			pos3D.Delete(i); int3D.erase(int3D.begin() + i); tempPredictions.Delete(i);	// delete the particle and
																						// if the particles' translation b/w frames is within the threshold of largestParticleShift
//			if (Distance(activeLongTracks[i].Last(), activeLongTracks[i].Penultimate()) <= shiftThreshold && Distance(activeLongTracks[i].Penultimate(), activeLongTracks[i].Antepenultimate()) <= shiftThreshold)
			if (activeLongTracks[i].Length() >= 7)
				exitTracks.push_back(activeLongTracks[i]);								// shift the corresponding activeLongTrack to exitTracks (or)
//			else {
//				inactiveTracks.push_back(activeLongTracks[i]);							// shift the corresponding activeLongTrack to inactiveTracks
//				a_is++;
//			}
			
			activeLongTracks.erase(activeLongTracks.begin() + i); s_al++;
			ignoreCam[i] = 100;
		}
		// TODO
		//else if (belowIntensity + leftcams >= 2) {					//  if the particle disappears on at least 2 cams (out of bounds or below intensity thresh)
		//	pos3D.Delete(i); int3D.erase(int3D.begin() + i); tempPredictions.Delete(i);	// delete the particle and
		//																				// if the track is longer than 12 frames and particles' translation b/w frames is within the threshold of largestParticleShift
		//	if (activeLongTracks[i].Length() >= 12 && Distance(activeLongTracks[i].Last(), activeLongTracks[i].Penultimate()) <= shiftThreshold) {
		//		inactiveLongTracks.push_back(activeLongTracks[i]);	s3++;				// shift the corresponding activeLongTrack to inactiveLongTracks (or)
		//	}
		//	else {
		//		inactiveTracks.push_back(activeLongTracks[i]);							// shift the corresponding activeLongTrack to inactiveTracks
		//		a3++;
		//	}

		//	activeLongTracks.erase(activeLongTracks.begin() + i); s2++;
		//	ignoreCam[i] = 100;
		//}
		else
			i++;
	}

	return ignoreCam;
}

// performing IPR on the residual images after removing tracked particles
Frame STB::IPRonResidual(Calibration& calib, Tiff2DFinder& t, deque<int**>& pixels_orig, deque<int**>& pixels_reproj, deque<int**>& pixels_res, Frame& estimates) {
	double mindist_3D = _ipr.mindist_3D;														// performing triangulation / IPR on the remaining particles in residual images
	double mindist_2D = _ipr.mindist_2D;
//	double xlim = ((double)(Npixw - 1)) / 2, ylim = ((double)(Npixh - 1)) / 2;
//	double intensityThresh = 0.8*_ipr.Get_threshold();//0.25*_ipr.Get_threshold();
	deque<Position> candidates;
	deque<int> camNums;
	for (int i = 0; i < ncams; i++)
		camNums.push_back(i);

	deque<Camera> cam_param = calib.Get_cam();

	double pix_length = (cam_param[0].Get_wpix() + cam_param[0].Get_hpix()) / 2; // this will be used as the shift increment

	for (int outerloop = 0; outerloop < _ipr.it_outerloop; outerloop++) {						// identify the particle candidates from residual images 

//		calib.Set_min2D(pow(1.1, outerloop)*mindist_2D);
		calib.Set_min2D(mindist_2D + pix_length / 4 * outerloop);

		// running IPR on the residual images
		Frame temp = _ipr.IPRLoop(calib, OTFcalib, camNums, ALL_CAMS, t.Get_colors(), pixels_orig, pixels_reproj, pixels_res, outerloop);
		candidates.insert(candidates.end(), temp.begin(), temp.end());
		cout << "\t# of particles detected in outerloop" << outerloop << ": " << temp.NumParticles() << endl;
		cout << "\tTotal particles (" << candidates.size() << ")" << endl;
	}

	if (_ipr.reducedCams) {																		// ipr with reduced cams

		Calibration calibReduced(_ipr.calibfile, REDUCED_CAMS, _ipr.mindist_2Dr, _ipr.mindist_3Dr, ncams);		// new calibration file for reduced cameras					
		for (int outerloop = 0; outerloop < _ipr.it_reducedCam; outerloop++) {

//			calibReduced.Set_min2D(pow(1.1, outerloop)*_ipr.mindist_2Dr);
			calibReduced.Set_min2D(_ipr.mindist_2Dr + pix_length / 4 * outerloop);

																								// running IPR by ignoring cameras one by one
			for (int ignoreCam = 0; ignoreCam < ncams; ignoreCam++) {
				Frame temp = _ipr.IPRLoop(calibReduced, OTFcalib, camNums, ignoreCam, t.Get_colors(), pixels_orig, pixels_reproj, pixels_res, outerloop);
//				deque<Position> temp_pos;
//				temp_pos.insert(temp_pos.end(), temp.begin(), temp.end());
//				_ipr.SaveParticlePositions(temp_pos, "/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/ParticlePositions/totalparticles.txt");
				//if (ignoreCam == 0)																// if the projection of a particle candidate on ignored cam is within its bounds and finds no particle there,
				//	for (int i = temp.NumParticles()-1; i >= 0 ; i--) 							// then delete the particle candidate						
				//		if (abs(temp[i].X1() - xlim) <= xlim && abs(temp[i].Y1() - ylim) <= ylim && pixels_orig[0][(int)round(temp[i].Y1())][(int)round(temp[i].X1())] < intensityThresh)
				//			temp.Delete(i);

				//if (ignoreCam == 1)																
				//	for (int i = temp.NumParticles() - 1; i >= 0; i--) 												
				//		if (abs(temp[i].X2() - xlim) <= xlim && abs(temp[i].Y2() - ylim) <= ylim && pixels_orig[1][(int)round(temp[i].Y2())][(int)round(temp[i].X2())] < intensityThresh)
				//			temp.Delete(i);

				//if (ignoreCam == 2)																
				//	for (int i = temp.NumParticles() - 1; i >= 0; i--) 												
				//		if (abs(temp[i].X3() - xlim) < xlim && abs(temp[i].Y3() - ylim) < ylim && pixels_orig[2][(int)round(temp[i].Y3())][(int)round(temp[i].X3())] < intensityThresh)
				//			temp.Delete(i);

				//if (ignoreCam == 3)															
				//	for (int i = temp.NumParticles() - 1; i >= 0; i--) 												
				//		if (abs(temp[i].X4() - xlim) < xlim && abs(temp[i].Y4() - ylim) < ylim && pixels_orig[3][(int)round(temp[i].Y4())][(int)round(temp[i].X4())] < intensityThresh)
				//			temp.Delete(i);
					
				candidates.insert(candidates.end(), temp.begin(), temp.end());
				cout << "\t# of particles detected in outerloop" << outerloop << ", ignoring cam" << ignoreCam << ": " << temp.NumParticles() << endl;
			}
			cout << "\tTotal particles (" << candidates.size() << ")" << endl;
		}

	}

	for (int i = 0; i < candidates.size(); i++) {												// removing a candidate if it's within 1 pixel of tracked particles
		double X1 = candidates[i].X1(), X2 = candidates[i].X2(), X3 = candidates[i].X3(), X4 = candidates[i].X4();
		double Y1 = candidates[i].Y1(), Y2 = candidates[i].Y2(), Y3 = candidates[i].Y3(), Y4 = candidates[i].Y4();

		for (int j = 0; j < estimates.NumParticles(); j++) {
			int num_overlap = 0;
			int dist_thred = 0.5;// _ipr.Get_psize() / 2;
			if (fabs(X1 - estimates[j].X1()) <= dist_thred && fabs(Y1 - estimates[j].Y1()) <= dist_thred) num_overlap ++;
			if (fabs(X2 - estimates[j].X2()) <= dist_thred && fabs(Y2 - estimates[j].Y2()) <= dist_thred) num_overlap ++;
			if (fabs(X3 - estimates[j].X3()) <= dist_thred && fabs(Y3 - estimates[j].Y3()) <= dist_thred) num_overlap ++;
			if (fabs(X4 - estimates[j].X4()) <= dist_thred && fabs(Y4 - estimates[j].Y4()) <= dist_thred) num_overlap ++;
//			if ((abs(X1 - estimates[j].X1()) < 1 && abs(Y1 - estimates[j].Y1()) < 1) &&
//				(abs(X2 - estimates[j].X2()) < 1 && abs(Y2 - estimates[j].Y2()) < 1) &&
//				(abs(X3 - estimates[j].X3()) < 1 && abs(Y3 - estimates[j].Y3()) < 1) &&
//				(abs(X4 - estimates[j].X4()) < 1 && abs(Y4 - estimates[j].Y4()) < 1))
			if (num_overlap >= 4) {
				candidates.erase(candidates.begin() + i);
				break;
			}
		}
	}
	cout << "\tTotal particles after removing close particles (" << candidates.size() << ")" << endl;



	return Frame(candidates);
}

// linking short tracks with particle candidates in residual images
void STB::MakeShortLinkResidual(int nextFrame, Frame& candidates, deque<Track>::iterator& tr, int iterations, bool* erase, bool* candidate_used) {

	pair<Frame::const_iterator, float> cost;

	for (int it = 0; it < iterations; it++) {													// iteratively trying to find a link for the short track from particle candidates
		double rsqr = pow(pow(1.1, it) * 3 * avgIPDist, 2);
		double shift = pow(1.1, it) * largestShift;
		deque<double> dist;
		double totalDist = 0;
		deque<Position> disp;
		Position vel(0, 0, 0);																	// calculating the predictive vel. field as an avg. of particle velocities from neighbouring tracks
		double d;

//#pragma omp parallel //shared(dist, disp)
//		{
//#pragma omp for
		for (int j = 0; j < activeLongTracks.size(); ++j) {										// identifying the neighbouring tracks (using 3*avg interparticle dist.) and getting their particle velocities
			double dsqr = Distance(tr->Last(), activeLongTracks[j].Penultimate());


			if (dsqr < rsqr && dsqr > 0) {
				d = pow(dsqr, 0.5);
				totalDist = totalDist + d;
//#pragma omp critical//(dist)
//				{
				dist.push_back(d);
				disp.push_back(activeLongTracks[j].Last() - activeLongTracks[j].Penultimate());
//				}
//			}
		}
		}

		if (dist.size() > 0) {															// if it finds neighbouring tracks
			deque<double> anti_dist;
			double max_dist = dist[0];
			for (int j = 0; j < dist.size(); j++) {
				if (max_dist < dist[j])
					max_dist = dist[j];
			}
			for (int j = 0; j < dist.size(); j++) {
				anti_dist.push_back(max_dist / dist[j]);
			}
			double total_anti_dist;
			for (int j = 0; j < dist.size(); j++) {
				total_anti_dist = total_anti_dist + anti_dist[j];
			}
			for (int j = 0; j < dist.size(); j++) {												// perform Gaussian wt. avg. to get velocity field
				double weight = (dist[j] / totalDist);
//				double weight = (anti_dist[j] / total_anti_dist);
				vel.Set_X(vel.X() + weight*disp[j].X());
				vel.Set_Y(vel.Y() + weight*disp[j].Y());
				vel.Set_Z(vel.Z() + weight*disp[j].Z());
			}

			Position estimate = tr->Last() + vel;
//			if (fabs(tr->Last().X() - 2.2104) < 0.01 ) {
//							int catch_label = 1;
//						}																					// finding a link for this short track using the velocity field
			cost = ComputeCost(candidates, candidates, searchRadiusSTB, estimate, vel, vel, true, candidate_used);

			if (cost.second == UNLINKED) {
				Position estimate = tr->Last();														// using nearest neighbour to make a link
				cost = ComputeCost(candidates, candidates, shift, estimate, vel, vel, true, candidate_used);
			}

		}

		else {																			// if no neighbouring tracks are identified
			Position estimate = tr->Last();														// using nearest neighbour to make a link
			cost = ComputeCost(candidates, candidates, shift, estimate, vel, vel, true, candidate_used);
		}

		if (cost.second != UNLINKED) {															// if linked with a candidate
			tr->AddNext(*cost.first, nextFrame);												// add the candidate to the track, delete it from list of untracked candidates and break the loop
//			candidate_used[cost.first.where()] = true;
//			candidates.Delete(cost.first.where());      // Shiyong Tan: Label them in ComputeCost instead of deleting them for the convenience of parallelization.
//			++tr;
			break;
		}
	}

	if (cost.second == UNLINKED) {																// if no link is found for the short track
		if (tr->Length() > 3) {																	// and if the track has at least 4 particles
//			inactiveTracks.push_back(*tr);														// add it to inactive tracks
			a_is++;
		}
		if (tr->Length() == 1)
			s_as1++;
		else if (tr->Length() == 2)
			s_as2++;
		else if (tr->Length() == 3)
			s_as3++;

//		tr = activeShortTracks.erase(tr);														// then delete from activeShortTracks
		*erase= true;  // Shiyong Tan: Label them instead of deleting them for the convenience of parallelization.
		
	}
}

// linking short tracks with particle candidates in residual images
void STB::MakeShortLinkResidual(int nextFrame, Frame& candidates, deque<Track>::iterator& tr, int iterations) {

	pair<Frame::const_iterator, float> cost;

	for (int it = 0; it < iterations; it++) {													// iteratively trying to find a link for the short track from particle candidates
		double rsqr = pow(pow(1.1, it) * 3 * avgIPDist, 2);
		double shift = pow(1.1, it) * largestShift;
		deque<double> dist;
		double totalDist = 0;
		deque<Position> disp;
		Position vel(0, 0, 0);																	// calculating the predictive vel. field as an avg. of particle velocities from neighbouring tracks
		double d;

//#pragma omp parallel //shared(dist, disp)
//		{
//#pragma omp for
		for (int j = 0; j < activeLongTracks.size(); ++j) {										// identifying the neighbouring tracks (using 3*avg interparticle dist.) and getting their particle velocities
			double dsqr = Distance(tr->Last(), activeLongTracks[j].Penultimate());
			if (dsqr < rsqr && dsqr > 0) {
				d = pow(dsqr, 0.5);
				totalDist = totalDist + d;
//#pragma omp critical//(dist)
//				{
				dist.push_back(d);
				disp.push_back(activeLongTracks[j].Last() - activeLongTracks[j].Penultimate());
//				}
//			}
		}
		}

		if (dist.size() > 0) {															// if it finds neighbouring tracks
			for (int j = 0; j < dist.size(); j++) {												// perform Gaussian wt. avg. to get velocity field
				double weight = (dist[j] / totalDist);
				vel.Set_X(vel.X() + weight*disp[j].X());
				vel.Set_Y(vel.Y() + weight*disp[j].Y());
				vel.Set_Z(vel.Z() + weight*disp[j].Z());
			}

			Position estimate = tr->Last() + vel;
																								// finding a link for this short track using the velocity field
			cost = ComputeCost(candidates, candidates, searchRadiusSTB, estimate, vel, vel, true);
		}

		else {																			// if no neighbouring tracks are identified
			Position estimate = tr->Last();														// using nearest neighbour to make a link
			cost = ComputeCost(candidates, candidates, shift, estimate, vel, vel, true);
		}

		if (cost.second != UNLINKED) {															// if linked with a candidate
			tr->AddNext(*cost.first, nextFrame);												// add the candidate to the track, delete it from list of untracked candidates and break the loop
			candidates.Delete(cost.first.where());
			++tr;
			break;
		}
	}

	if (cost.second == UNLINKED) {																// if no link is found for the short track
		if (tr->Length() > 3) {																	// and if the track has at least 4 particles
//			inactiveTracks.push_back(*tr);														// add it to inactive tracks
			a_is++;
		}
		if (tr->Length() == 1)
			s_as1++;
		else if (tr->Length() == 2)
			s_as2++;
		else if (tr->Length() == 3)
			s_as3++;

		tr = activeShortTracks.erase(tr);														// then delete from activeShortTracks

	}
}

bool STB::CheckVelocity(Track& tr) {
	int num_track_to_average = 25;
	double search_radius = designated_search_radius; // TODO: make it as a configure parameter
//	double threshold_scale = 5;
	if (particle_search_radius > 0) search_radius = particle_search_radius;
	double threshold_scale = 5 + 10 * search_radius / designated_search_radius;

//	vector<double> dist;
//	vector<int> neighbor_track;
	double dist[num_track_to_average];
	int neighbor_track[num_track_to_average];

	Position last = tr.Last();
	double x_l = last.X() - search_radius, x_u = last.X() + search_radius;
	double y_l = last.Y() - search_radius, y_u = last.Y() + search_radius;
	double z_l = last.Z() - search_radius, z_u = last.Z() + search_radius;

	int num_track = 0;

#pragma omp parallel shared(num_track) //num_threads(8)
						{
#pragma omp for
	for (int i = 0; i < activeLongTracks.size(); ++i) {
		double x_test = activeLongTracks[i].Last().X();

//		double d = pow(Distance(tr->Last(), activeLongTracks[i].Last()), 0.5);
//		if (fabs(tr.Last().X() - activeLongTracks[i].Last().X()) < search_radius)
//			if (fabs(tr.Last().Y() - activeLongTracks[i].Last().Y()) < search_radius)
//				if (fabs(tr.Last().Z() - activeLongTracks[i].Last().Z()) < search_radius) {
		if (x_test > x_l && x_test < x_u) {
			double y_test = activeLongTracks[i].Last().Y();
				if (y_test > y_l && y_test < y_u) {
					double z_test = activeLongTracks[i].Last().Z();
					if (z_test > z_l && z_test < z_u ) {
						double d = pow(Distance(tr.Last(), activeLongTracks[i].Last()), 0.5);
	//		if (d < search_radius) {
	//			if (neighbor_track.size() <= 50) {
	#pragma omp critical(num_track)
								{
				if (num_track < num_track_to_average) {
					dist[num_track] = d;
					neighbor_track[num_track] = i;
					++ num_track;
	//				dist.push_back(d);
	//				neighbor_track.push_back(i);
				} else { // if the number is larger than 50, then keep the closest 50
					// get the maximun dist
					double dist_max = dist[0];
					int index_max = 0;
	//				for (int j = 1; j < neighbor_track.size(); ++j) {
					for (int j = 1; j < num_track_to_average; ++j) {
						if (dist[j] > dist_max) {
							dist_max = dist[j];
							index_max = j;
						}
					}
					if (d < dist_max) {
	//					dist.erase(dist.begin() + index_max);
	//					dist.push_back(d);
	//					neighbor_track.erase(neighbor_track.begin() + index_max);
	//					neighbor_track.push_back(i);
						dist[index_max] = d;
						neighbor_track[index_max] = i;
					}
				}
								}
					}
				}
		}
	}
						}

//	if (neighbor_track.size() < 10) return true; // if less than 10 track is found, then return true
	if (num_track < 0) return true;

	if (num_track == num_track_to_average) {
		double mean_radius = 0;
		for (int i = 1; i < num_track_to_average; ++i) {
			mean_radius += dist[i];
		}
		mean_radius /= num_track_to_average;
		particle_search_radius = mean_radius;
	}

	// get the minimum radius
//	double min_radius =  dist[0];
//	for (int i = 1; i < neighbor_track.size(); ++i) {
//		if (dist[i] < min_radius) {
//			min_radius = dist[i];
//		}
//	}
//	particle_dist = min_radius;

	vector<double> velocity;
//	for (int i = 0; i < neighbor_track.size(); ++i) {
	for (int i = 0; i < num_track; ++i) {
		velocity.push_back(pow(Distance(activeLongTracks[neighbor_track[i]].Last(), activeLongTracks[neighbor_track[i]].Penultimate()), 0.5));
	}

	double mean_velocity = 0;
	for (int i = 0; i < velocity.size(); ++i) {
		mean_velocity = mean_velocity + velocity[i];
	}
	mean_velocity = mean_velocity / velocity.size();

	double diff_velocity = 0;
	for (int i = 0; i < velocity.size(); ++i) {
		diff_velocity = diff_velocity + pow((velocity[i] - mean_velocity), 2);
	}
	diff_velocity = pow(diff_velocity / velocity.size(), 0.5);

	double velocity_self = pow(Distance(tr.Last(), tr.Penultimate()), 0.5);

	if (fabs(velocity_self - mean_velocity) > threshold_scale * diff_velocity) return false;

	return true;
}

bool STB::CheckAcceleration(deque<Track>::iterator& tr) {
	if (!enable_acc_check) return true; // if not enable acceleration check then return true.

	int len = tr->Length();
	Position point0 = tr->GetPos(len - 1);
	Position point1 = tr->GetPos(len - 2);
	Position point2 = tr->GetPos(len - 3);
	Position point3 = tr->GetPos(len - 4);

	Position vel_1 = point0 - point1;
	Position vel_2 = point1 - point2;
	Position vel_3 = point2 - point3;

	double acc_1 = pow(Distance(vel_1, vel_2), 0.5);
	double acc_2 = pow(Distance(vel_2, vel_3), 0.5);

	if (fabs(acc_1 - acc_2) > acc_diff_thred) return false;

	return true;
}

void STB::GetAccThred() {
	vector<double> acc_diff_vec;
	for (int i = 0; i < activeLongTracks.size(); ++i) {
		int len = activeLongTracks[i].Length();
		if (len < 20) continue; // only tracks of which the length is longer than 20 are considered as good tracks
		int num_frame_to_average = 5; // average the acceleration difference for the last 10 frames
		double acc_diff = 0;
		for (int j = len - num_frame_to_average; j < len; ++j) {
			Position point0 = activeLongTracks[i].GetPos(j - 1);
			Position point1 = activeLongTracks[i].GetPos(j - 2);
			Position point2 = activeLongTracks[i].GetPos(j - 3);
			Position point3 = activeLongTracks[i].GetPos(j - 4);

			Position vel_1 = point0 - point1;
			Position vel_2 = point1 - point2;
			Position vel_3 = point2 - point3;

			double acc_1 = pow(Distance(vel_1, vel_2), 0.5);
			double acc_2 = pow(Distance(vel_2, vel_3), 0.5);

			acc_diff = acc_diff + fabs(acc_1 - acc_2);
		}
		acc_diff = acc_diff / num_frame_to_average;
		acc_diff_vec.push_back(acc_diff);
	}
	if (acc_diff_vec.size() > 100) { // only calculate the threshold when data is sufficient
		double acc_diff_ave = 0;
		for (int i = 0; i < acc_diff_vec.size(); ++i) {
			acc_diff_ave = acc_diff_ave + acc_diff_vec[i];
		}
		acc_diff_ave = acc_diff_ave / acc_diff_vec.size();

		double acc_diff_std = 0;
		for (int i = 0; i < acc_diff_vec.size(); ++i) {
			acc_diff_std = acc_diff_std + pow((acc_diff_vec[i] - acc_diff_ave), 2);
		}
		acc_diff_std = pow(acc_diff_std / acc_diff_vec.size(), 0.5);

		acc_diff_thred = acc_diff_ave + 10 * acc_diff_std;

		enable_acc_check = true;
	} else {
		enable_acc_check = false;
	}
}

void SimpleLinearRegression(double* y, double* coeff) {
	int x[4];
	for (int i = 0; i < 4; ++i) x[i] = i + 1;
    double xsum = 0,x2sum = 0,ysum = 0,xysum = 0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
    for (int i=0; i < 4; ++i) {
        xsum = xsum + x[i];                        //calculate sigma(xi)
        ysum = ysum + y[i];                        //calculate sigma(yi)
        x2sum = x2sum + pow(x[i], 2);                //calculate sigma(x^2i)
        xysum = xysum + x[i] * y[i];                    //calculate sigma(xi*yi)
    }
    coeff[0] = (4 * xysum - xsum * ysum) / (4 * x2sum - xsum * xsum);            //calculate slope
    coeff[1] = (x2sum * ysum - xsum * xysum) / (x2sum * 4 - xsum * xsum);            //calculate intercept

}

//bool STB::CheckLinearFit(deque<Track>::iterator& tr) {
bool STB::CheckLinearFit(Track& tr) {
	int len = tr.Length();

//	Position point0 = tr.GetPos(len - 4);
//	Position point1 = tr.GetPos(len - 3);
//	Position point2 = tr.GetPos(len - 2);
//	Position point3 = tr.GetPos(len - 1);

	Position point[4];
	for (int i = 0; i < 4; ++ i) {
		point[i] = tr.GetPos(len - 4 + i);
	}

	double y[4];
	y[0] = point[0].X(); y[1] = point[1].X(); y[2] = point[2].X(); y[3] = point[3].X();
	double coeffx[2];
	SimpleLinearRegression(y, coeffx);
	double x_fit = coeffx[0] * 4 + coeffx[1];

	y[0] = point[0].Y(); y[1] = point[1].Y(); y[2] = point[2].Y(); y[3] = point[3].Y();
	double coeffy[2];
	SimpleLinearRegression(y, coeffy);
	double y_fit = coeffy[0] * 4 + coeffy[1];

	y[0] = point[0].Z(); y[1] = point[1].Z(); y[2] = point[2].Z(); y[3] = point[3].Z();
	double coeffz[2];
	SimpleLinearRegression(y, coeffz);
	double z_fit = coeffz[0] * 4 + coeffz[1];

	Position point_fit(x_fit, y_fit, z_fit);

	double del_pos[4];
	del_pos[3] = pow(Distance(point[3], point_fit), 0.5);

	double error_thred = .05;
//	if (linear_fit_error[99] == 1) {
//		error_thred = linear_fit_error[0];
//	}

	if (del_pos[3] > error_thred) return false;

	if (del_pos[3] > error_thred * 0.6) {  //For those tracks that are at the edge of being deleted, check the whole track to see whether its error is large also
		double total_del = del_pos[3];
		for (int i = 0; i < 3; ++ i) {
			x_fit = coeffx[0] * i + coeffx[1];
			y_fit = coeffy[0] * i + coeffy[1];
			z_fit = coeffz[0] * i + coeffz[1];
			Position point_fit(x_fit, y_fit, z_fit);
			del_pos[i] = pow(Distance(point[i], point_fit), 0.5);
			if (del_pos[i] > error_thred * 0.6) return false;
			total_del += del_pos[i];
		}
		if (total_del / 4 > error_thred * 0.6) return false;
	}

//	if (linear_fit_error[99] == 0) {
//		int num_sample =0;
//		for (int i = 0; i < 99; ++i) {
//			if (linear_fit_error[i] == 0) {
//				linear_fit_error[i] = del_pos;
//				break;
//			}
//			++ num_sample;
//		}
//		if (num_sample == 99) {
//			double fit_error = 0;
//			for (int i = 0; i < 99; ++i) {
//				fit_error += linear_fit_error[i];
//			}
//			fit_error /= 99;
//			double fit_std = 0;
//			for (int i = 0; i < 99; ++i) {
//				fit_std += pow(fit_error - linear_fit_error[i], 2);
//			}
//			fit_std = pow(fit_std / 99, .5);
//			linear_fit_error[0] = fit_error + fit_std * 5;
//			linear_fit_error[99] = 1;
//		}
//	}

	return true;

}



//########################### MAT / TXT FILES ###################################

void STB::MatTracksSave(string address, string s, bool is_back_STB) {
	// Saving tracks for Matlab
	string X1 = "ActiveLongTracks" + s, X2 = "ActiveShortTracks" + s, X3 = "InactiveTracks" + s, X4 = "ExitTracks" + s, X5 = "InactiveLongTracks" + s;
	
/*
 * Modified By Shiyong Tan, 8/12/18
 * The previous saving method has a lot of zeros points, which make the file very large
 * Thus, change the saving format.
 * Start:
 */
//	MatfileSave(activeLongTracks, address + X1, X1, lastFrame);
//	MatfileSave(activeShortTracks, address + X2, X2, lastFrame);
//	MatfileSave(inactiveTracks, address + X3, X3, lastFrame);
//	MatfileSave(exitTracks, address + X4, X4, lastFrame);
//	MatfileSave(inactiveLongTracks, address + X5, X5, lastFrame);
	
	SaveTrackToTXT(activeLongTracks, address + X1);
	SaveTrackToTXT(activeShortTracks, address + X2);
//	SaveTrackToTXT(inactiveTracks, address + X3);  // No longer save inactiveTracks which are useless
	SaveTrackToTXT(exitTracks, address + X4);
	SaveTrackToTXT(inactiveLongTracks, address + X5);
	if (is_back_STB) {
		SaveTrackToTXT(bufferTracks, address + "BufferTracks" + s);
	}
// End
}

void STB::LoadAllTracks(string address, string frame_number, bool is_back_STB) {
	string s = frame_number + ".txt";
	string X1 = "ActiveLongTracks" + s, X2 = "ActiveShortTracks" + s, X3 = "InactiveTracks" + s, X4 = "ExitTracks" + s, X5 = "InactiveLongTracks" + s, X6 = "BufferTracks" + s;
/*
 * Modified by Shiyong Tan, 8/12/18
 * Adapt to the change of the saving format.
 * Start:
 */
//	Load_Tracks(address + X1, ActiveLong);
//	Load_Tracks(address + X2, ActiveShort);
//	Load_Tracks(address + X3, Inactive);
//	Load_Tracks(address + X4, Exit);
//	Load_Tracks(address + X5, InactiveLong);
//
	LoadTrackFromTXT(address + X1, ActiveLong);
	LoadTrackFromTXT(address + X2, ActiveShort);
//	LoadTrackFromTXT(address + X3, Inactive);  // no longer load inactive long tracks, which are useless
	if (error == NONE) {
		LoadTrackFromTXT(address + X4, Exit);
		LoadTrackFromTXT(address + X5, InactiveLong);
		if (is_back_STB) LoadTrackFromTXT(address + X6, Buffer);
		error = NONE;
	} else {
		LoadTrackFromTXT(address + X4, Exit);
		LoadTrackFromTXT(address + X5, InactiveLong);
		if (is_back_STB) LoadTrackFromTXT(address + X6, Buffer);
	}

// End
}


void STB::MatfileSave(deque<Track> tracks, string address, string name, int size) {
/*
 * Modified by Shiyong Tan, 2/8/18
 * Discard using matio, use Data_IO instead
 * Start:
 */
//	// Create a .mat file with pos3D
//	size_t sizeofpos3D = tracks.size();
//	double* tempX = new double[size];
//	double* tempY = new double[size];
//	double* tempZ = new double[size];
//	mat_t    *matfp;
//	matvar_t *cell_arrayX, *cell_elementX, *cell_arrayY, *cell_elementY, *cell_arrayZ, *cell_elementZ;
//	size_t dims[2] = { sizeofpos3D, 1 };
//	string mat_nameX = name + "X";
//	string mat_nameY = name + "Y";
//	string mat_nameZ = name + "Z";
//	matfp = Mat_CreateVer((address +".mat").c_str(), NULL, MAT_FT_DEFAULT);
//
//	switch (NULL == matfp) {
//		fprintf(stderr, "Error creating MAT file \"Tracking\"!\n");
//		break;
//	}
//	cell_arrayX = Mat_VarCreate(mat_nameX.c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//	cell_arrayY = Mat_VarCreate(mat_nameY.c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//	cell_arrayZ = Mat_VarCreate(mat_nameZ.c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//	if (NULL == cell_arrayX || NULL == cell_arrayY || NULL == cell_arrayZ) {
//		fprintf(stderr, "Error creating variable for 'Tracking'\n");
//	}
//	else {
//		for (int i = 0; i < sizeofpos3D; i++) {
//			dims[0] = 1;
//			dims[1] = size;
//			int time = 0;
//			for (int k = 0; k < size; k++) {
//				if (time < tracks[i].Length() && k == tracks[i].GetTime(time) - 1) {
//					tempX[k] = tracks[i][time].X();
//					tempY[k] = tracks[i][time].Y();
//					tempZ[k] = tracks[i][time].Z();
//					time++;
//				}
//				else {
//					tempX[k] = 0;
//					tempY[k] = 0;
//					tempZ[k] = 0;
//				}
//			}
//
//			cell_elementX = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tempX, 0);
//			cell_elementY = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tempY, 0);
//			cell_elementZ = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tempZ, 0);
//			switch (NULL == cell_elementX || NULL == cell_elementY || NULL == cell_elementZ) {
//				fprintf(stderr, "Error creating cell element variable\n");
//				Mat_VarFree(cell_arrayX); //Mat_VarFree(cell_arrayY); Mat_VarFree(cell_arrayZ);
//				Mat_Close(matfp);
//				break;
//			}
//			Mat_VarSetCell(cell_arrayX, i, cell_elementX);
//			Mat_VarSetCell(cell_arrayY, i, cell_elementY);
//			Mat_VarSetCell(cell_arrayZ, i, cell_elementZ);
//		}
//	}
//
//	Mat_VarWrite(matfp, cell_arrayX, MAT_COMPRESSION_NONE);
//	Mat_VarWrite(matfp, cell_arrayY, MAT_COMPRESSION_NONE);
//	Mat_VarWrite(matfp, cell_arrayZ, MAT_COMPRESSION_NONE);
//	Mat_VarFree(cell_arrayX); Mat_VarFree(cell_arrayY); Mat_VarFree(cell_arrayZ);
//	Mat_Close(matfp);
//	delete[] tempX, tempY, tempZ;
	// TODO: to check whether it works.

	// convert track into 3D matrix the time for each track is the same.
	size_t sizeofpos3D = tracks.size(); // the number of particle
	/*
	 * Modified by Shiyong Tan
	 * date: 5/22/18
	 * !!!!IMPORTANT!!!!
	 *  it is a very difficult bug
	 * double track_data will be saved in stack. Its size is unknown.
	 * when its size is large, it may overwrite memory of other variable in the stack.
	 * It is illegal to use it. But don't know why the compiler let it pass.
	 * A way to avoid it is to put this data into heap instead of stack, so that it won't overwrite other variable.
	 * Start:
	 */
//	double track_data[sizeofpos3D][size][3];  // size is the total number of frames
	double* track_data = new double[sizeofpos3D * size * 3];
	//END

	for(int i = 0; i < sizeofpos3D; i++) {
		int absolute_starttime = 0;
		int absolute_endtime = 0;
		absolute_starttime = tracks.at(i).GetTime(0); // 0 is the starting time of the begining of the motion of a particle
										// absolute start time is the start time of a particle in the overall time reference for all particles.
		absolute_endtime = tracks.at(i).GetTime(tracks.at(i).Length() - 1);
										//tracks[i].Length() is the end of the motion of a particle
										// absolute end time is the end time of a particle in the overall time reference for all particle
		int time = 0;
		for(int j = 0; j < size; j++) {
			if (absolute_starttime <= j + 1 && j + 1 <= absolute_endtime) {
				track_data[i * size * 3 + j * 3 + 0] = tracks.at(i)[time].X();
				track_data[i * size * 3 + j * 3 + 1] = tracks.at(i)[time].Y();
				track_data[i * size * 3 + j * 3 + 2] = tracks.at(i)[time].Z();
				time ++;
			} else { // the frame the particle doesn't show up is set as 0.
				track_data[i * size * 3 + j * 3 + 0] = 0;
				track_data[i * size * 3 + j * 3 + 1] = 0;
				track_data[i * size * 3 + j * 3 + 2] = 0;
			}
		}
	}

	double dimension_info[3]; // a vector to save dimension info, this vector will tell how to read data
	dimension_info[0] = sizeofpos3D; // number of elements in 1st dimension
	dimension_info[1] = size; // number of elements in 2nd dimension
	dimension_info[2] = 3; // number of elements in 3rd dimension

	NumDataIO<double> data_io;
	data_io.SetFilePath(address + ".txt");

	// to save dimension info
	data_io.SetTotalNumber(3);
	data_io.SetNumPrecsion(12);
	data_io.WriteData((double*) dimension_info);

	// to append data
	data_io.SaveMode(1); //to set the save mode as append
	data_io.SetTotalNumber(sizeofpos3D * size * 3);
	data_io.WriteData((double*) track_data);

	delete[] track_data;

	// End
}

void STB::SaveTrackToTXT(deque<Track> tracks, string address) {
	size_t num_track = tracks.size(); // number of tracks
	unsigned int total_len = 0; // Total length of all the tracks
	for (unsigned int i  = 0; i < num_track; ++ i) {
		total_len = total_len + tracks.at(i).Length();
	}
	double* track_data = new double[total_len * 5];

	unsigned int index = 0;
	for (unsigned int i = 0; i <  num_track; ++ i) {
		unsigned int len = tracks.at(i).Length();
		unsigned int start_time = tracks.at(i).GetTime(0);
		for (unsigned int j = 0; j < len; ++ j) {
			track_data[index] = i; // the track NO.
			++ index;
			track_data[index] = j + start_time; // the current frame
			++ index;
			track_data[index] = tracks.at(i)[j].X();
			++ index;
			track_data[index] = tracks.at(i)[j].Y();
			++ index;
			track_data[index] = tracks.at(i)[j].Z();
			++ index;
		}
	}

	NumDataIO<double> data_io;
	data_io.SetFilePath(address + ".txt");
	data_io.SetTotalNumber(total_len * 5);
	data_io.WriteData((double*) track_data);

	delete[] track_data;
}

//###################### TEMPORARY FUNTIONS FOR TESTING ###############################

// to load 3D positions without IPR
Frame STB::Load_3Dpoints(string path) {
	/*
	 * Modified by Shiyong Tan, 2/8/18
	 * Discard using matio, use Data_IO instead
	 * Start:
	 */
	Frame iprFrame;
//	//string file = "S:/Projects/Bubble/10.31.17/Run1/BubblesNParticlesLow_4000fps/ParticlesOnly/" + path + ".mat";
//	string file = tiffaddress + path + ".mat";
//	const char *fileName = file.c_str();
//	mat_t *mat = Mat_Open(fileName, MAT_ACC_RDONLY);
//
//	if (mat == NULL) {
//		cout << " error in reading the mat file" << endl;
//	}
//	string cam;
//	string name;
//
//	if (mat) {
//		//std::cout << "Open file to read\n\tmat == " << mat << "\n";
//
//		matvar_t *matVar = 0;
//		name = path;
//		matVar = Mat_VarRead(mat, (char*)name.c_str());
//
//		if (matVar) {
//			int rows;	int cols;
//			rows = matVar->dims[0]; cols = matVar->dims[1];
//			unsigned namesize = matVar->nbytes / matVar->data_size;
//			double *namedata = static_cast<double*>(matVar->data);
//			for (int i = 0; i < rows; i++) {
//				Position pos(namedata[i], namedata[rows + i], namedata[2*rows + i]);
//				iprFrame.Add(pos);
//			}
//		}
//	}
//	else {
//		cout << "Cannot open file\n";
//	}
//
//	Mat_Close(mat);
	// TODO: to check whether it works.
	string file = tiffaddress + path + ".txt";
	NumDataIO<double> data_io;
	data_io.SetFilePath(file);
	int total_num = data_io.GetTotalNumber();
	int rows = total_num / 3;
	double* points_array =  new double[rows * 3];
	data_io.ReadData(points_array);
	for (int i = 0; i < rows; i++) {
		Position pos(points_array[i * 3], points_array[i * 3 + 1], points_array[i * 3 + 2]);
		iprFrame.Add(pos);
	}
	return iprFrame;
	// END
}

// Loading the tracks from .mat file (For Testing)
void STB::Load_Tracks(string path, TrackType trackType) {
	/*
	 * Modified by Shiyong Tan, 2/8/18
	 * Discard using matio, use Data_IO instead
	 * Start:
	 */
//	cout << "Loading the tracks from .mat files" << endl;
//	string file = tiffaddress + path + ".mat";
//	const char *fileName = file.c_str();
//	mat_t *mat = Mat_Open(fileName, MAT_ACC_RDONLY);
//
//	if (mat == NULL) {
//		cout << " error in reading the mat file" << endl;
//	}
//
//	string namex, namey, namez;
//
//	if (mat) {
//		//std::cout << "Open file to read\n\tmat == " << mat << "\n";
//
//		matvar_t *matVarx = 0, *matVary = 0, *matVarz = 0;
//		namex = path + "X"; namey = path + "Y"; namez = path + "Z";
//		matVarx = Mat_VarRead(mat, (char*)namex.c_str());
//		matVary = Mat_VarRead(mat, (char*)namey.c_str());
//		matVarz = Mat_VarRead(mat, (char*)namez.c_str());
//
//		if (matVarx && matVary && matVarz) {
//			int rowsx = matVarx->dims[0], colsx = matVarx->dims[1];
//			double *namedatax = static_cast<double*>(matVarx->data);
//
//			int rowsy = matVary->dims[0], colsy = matVary->dims[1];
//			double *namedatay = static_cast<double*>(matVary->data);
//
//			int rowsz = matVarz->dims[0], colsz = matVarz->dims[1];
//			double *namedataz = static_cast<double*>(matVarz->data);
//
//			for (int i = 0; i < rowsx; i++) {
//				int time = 1;
//				Track track;
//				for (int j = 0; j < colsx; j++) {
//					int element = i + j*rowsx;
//					if (namedatax[element] == 0) {
//						time++; continue;
//					}
//
//					Position pos(namedatax[element], namedatay[element], namedataz[element]);
//					track.AddNext(pos, time);
//					time++;
//				}
//
//				switch (trackType)
//				{
//					case ActiveLong : activeLongTracks.push_back(track);
//										break;
//					case ActiveShort: activeShortTracks.push_back(track);
//										break;
//					case Inactive: inactiveTracks.push_back(track);
//										break;
//					case Exit: exitTracks.push_back(track);
//										break;
//					case InactiveLong: inactiveLongTracks.push_back(track);
//										break;
//				}
//			}
//		}
//	}
//	else {
//		cout << "Cannot open file\n";
//	}
//
//	Mat_Close(mat);

	NumDataIO<double> data_io;
	data_io.SetFilePath(path);

	//read the dimension info
	double dimension_info[3];
	data_io.SetTotalNumber(3);
	data_io.ReadData((double*) dimension_info);

	int num_particles = (int) dimension_info[0];
	int num_frames = (int) dimension_info[1];
	/*
		 * Modified by Shiyong Tan
		 * date: 5/22/18
		 * declare an unknown variable in stack should be forbidden
		 * Start:
	 */
//	double track_data[num_particles][num_frames][3];
	double* track_data = new double[num_particles * num_frames * 3];
	if (num_particles != 0)  {
		data_io.SetTotalNumber(num_particles * num_frames * 3);
			data_io.SetSkipDataNum(3); // skip the dimension info
			data_io.ReadData((double*) track_data);

			// convert 3D data into track
			for (int i = 0; i < num_particles; i++) {
				Track track;
				int time = 1;
				for (int j = 0; j < num_frames; j++) {
					if (track_data[i * num_frames * 3 + j * 3 + 0] == 0 &&
						track_data[i * num_frames * 3 + j * 3 + 1] == 0 &&
						track_data[i * num_frames * 3 + j * 3 + 2] == 0 ) { // if all of the data is 0, that means there is no track
						time++; continue;
					}
					Position pos(track_data[i * num_frames * 3 + j * 3 + 0],
							track_data[i * num_frames * 3 + j * 3 + 1], track_data[i * num_frames * 3 + j * 3 + 2]);
					track.AddNext(pos, time);
					time++;
				}
				switch (trackType)
				{
					case ActiveLong : activeLongTracks.push_back(track); break;
					case ActiveShort: activeShortTracks.push_back(track);break;
					case Inactive: inactiveTracks.push_back(track);break;
					case Exit: exitTracks.push_back(track);break;
					case InactiveLong: inactiveLongTracks.push_back(track);break;
				}
		}
	}
	delete[] track_data;
	//END
}

void STB::LoadTrackFromTXT(string path, TrackType trackType) {
	NumDataIO<double> data_io;
	data_io.SetFilePath(path);
	unsigned int total_num = data_io.GetTotalNumber();
	unsigned int total_len = total_num / 5;
	data_io.SetTotalNumber(total_num);
	double* track_data = new double[total_num];
	if (total_len != 0) {
		data_io.ReadData((double*) track_data);

		for (unsigned int i = 0; i < total_num; ) {
			unsigned int track_no = track_data[i];
			Track track;
			while (track_data[i] == track_no) { // for the same track
				++ i;
				int time = track_data[i];
				++ i;
				double X = track_data[i];
				++ i;
				double Y = track_data[i];
				++ i;
				double Z = track_data[i];
				++ i;
				Position pos(X, Y, Z);
				track.AddNext(pos, time);
			}
			// after getting the track
			switch (trackType)
			{
				case ActiveLong : activeLongTracks.push_back(track); break;
				case ActiveShort: activeShortTracks.push_back(track);break;
				case Inactive: inactiveTracks.push_back(track);break;
				case Exit: exitTracks.push_back(track);break;
				case InactiveLong: inactiveLongTracks.push_back(track);break;
				case Buffer: bufferTracks.push_back(track);break;
			}
		}
	}
}

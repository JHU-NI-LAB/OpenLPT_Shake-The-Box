#include <STB.h>
#include <BackSTB.h>
#include <chrono>
#include <cstdio>

#define ALL_CAMS 100
#define STBflag true
#define IPRflag false

using namespace std;

BackSTB::BackSTB(int firstFrame, int lastFrame, double threshold) : first(firstFrame), last(lastFrame), distThreshold(pow(threshold,2)) {}

void BackSTB::UpdateTracks(STB& s) {
	cout << "\tPerforming BackSTB" << endl;
	s.activeShortTracks.clear();
	auto start = std::chrono::system_clock::now();

	ncams = s.ncams;
	Npixh = s.Npixh;
	Npixw = s.Npixw;
	string address = s.tiffaddress + "Tracks/BackSTBTracks/";
	int frame_start = s.last;
	bool debug_load_success = false;
	if (debug_mode == SKIP_PREVIOUS_BACK_STB) {
		s.LoadAllTracks(address, to_string(debug_frame_number), 1);
		if (error == NO_FILE) {
			cout<<"No track files for the specific frame is found\n";
			error = NONE;
			s.LoadAllTracks(s.tiffaddress + "Tracks/ConvergedTracks/", to_string(s.last), 0);
		} else {
			frame_start = debug_frame_number + 4;
			debug_mode = NO_SKIP; //when execute the debug requirement, then label it as done.
			debug_load_success = true;
		}
	}
	// creating & initializing original, residual and reprojected image pixels
	for (int n = 0; n < s.ncams; n++) {
		pixels_orig.push_back(new int*[s.Npixh]);
		pixels_reproj.push_back(new int*[s.Npixh]);
		pixels_res.push_back(new int*[s.Npixh]);

		for (int i = 0; i < s.Npixh; i++) {
			pixels_orig[n][i] = new int[s.Npixw] {};
			pixels_reproj[n][i] = new int[s.Npixw] {};
			pixels_res[n][i] = new int[s.Npixw] {};
		}
	}

	Calibration calib(s._ipr.Get_calibfile(), ALL_CAMS, s._ipr.Get_mindist2D(), s._ipr.Get_mindist3D(), ncams);

	if (!debug_load_success) { // if the data is loaded from previous back STB , then it doesn't need to be preprocessed.
		// Move all tracks to active long new tracks
		for (deque<Track>::iterator tr = s.activeLongTracks.begin(); tr != s.activeLongTracks.end(); ) {
			s.bufferTracks.push_back(*tr);
			tr = s.activeLongTracks.erase(tr);
		}

		for (deque<Track>::iterator tr = s.exitTracks.begin(); tr != s.exitTracks.end(); ) {
			s.bufferTracks.push_back(*tr);
			tr = s.exitTracks.erase(tr);
		}

		for (deque<Track>::iterator tr = s.inactiveLongTracks.begin(); tr != s.inactiveLongTracks.end(); ) {
			s.bufferTracks.push_back(*tr);
			tr = s.inactiveLongTracks.erase(tr);
		}

		// trim the tracks
		cout<<"Trimming the tracks..." <<endl;
		for (deque<Track>::iterator tr = s.bufferTracks.begin(); tr != s.bufferTracks.end(); ) {
			if (tr->Length() > 15) {
				tr->DeleteFront(); tr->DeleteFront(); tr->DeleteFront();tr->DeleteFront();
				if (tr->GetTime(tr->Length() - 1) < s.last - 7) {
					tr->DeleteBack(); tr->DeleteBack();tr->DeleteBack();tr->DeleteBack();
				}
				++ tr;
			} else {
				tr = s.bufferTracks.erase(tr); // delete the tracks which are less than 15 particles
			}

		}
	}
	for (int currFrame = frame_start - 3; currFrame != first; currFrame--) {						// Linking the tracks back in time
//		clock_t start0, start1, start2;
//		start = clock();
			int prevFrame = currFrame - 1;
			if (prevFrame % 500 == 0) {
				cout<<"Loading tracks..."<<endl;
				// load inactive long tracks and exit tracks to buffer tracks
				// get the start of new tracks
//				deque<Track>::iterator tr_tr = s.bufferTracks.end(); // the address of the element after the end of the deque.
				int index_old_end = s.bufferTracks.size();
				s.LoadTrackFromTXT(s.tiffaddress + "Tracks/ConvergedTracks/InactiveLongTracks" + to_string(prevFrame) + ".txt", s.Buffer);
				s.LoadTrackFromTXT(s.tiffaddress + "Tracks/ConvergedTracks/ExitTracks" + to_string(prevFrame) + ".txt", s.Buffer);
				//  trim the tracks
	//			++ tr_tr;
				cout<<"Trimming the tracks..." <<endl;
				int index_new_end = s.bufferTracks.size();
				deque<Track>::iterator tr_tr = s.bufferTracks.end() - (index_new_end - index_old_end);
				for (; tr_tr != s.bufferTracks.end(); ) {
					if (tr_tr->Length() > 15) {
						tr_tr->DeleteFront(); tr_tr->DeleteFront(); tr_tr->DeleteFront();tr_tr->DeleteFront();
						tr_tr->DeleteBack(); tr_tr->DeleteBack();tr_tr->DeleteBack();tr_tr->DeleteBack();
						++ tr_tr;
					} else {
						tr_tr = s.bufferTracks.erase(tr_tr); // delete the tracks which are less than 15 particles
					}
				}
			}
			cout << "\tBackSTB at frame: " << prevFrame << endl;
			cout << "\t\tCheck Points: " << endl;

	//		start = clock();
			start =std::chrono::system_clock::now();
			cout << "\t\t\tLinking Long/Exit tracks: "<<endl;
			deque<Track> unlinkedTracks;
			Frame predictions; deque<double> intensity;
																								// for each track in activeLongTracks
	//		for (deque<Track>::iterator Ltr = s.activeLongTracks.begin(); Ltr != s.activeLongTracks.end(); ++Ltr)
	//			if (!(Ltr->Exists(prevFrame)) && Ltr->Exists(currFrame))					// if the track has no point in the prevFrame and has a point in the currFrame (unlinked track)
	//				Predictor(s, prevFrame, Ltr, unlinkedTracks, predictions, intensity);		// predict its position in prevFrame (using 4 time steps after prevFrame)
	//
	//																							// repeating the same for exit tracks
	//		for (deque<Track>::iterator Etr = s.exitTracks.begin(); Etr != s.exitTracks.end(); ++Etr)
	//			if (!(Etr->Exists(prevFrame)) && Etr->Exists(currFrame))
	//				Predictor(s, prevFrame, Etr, unlinkedTracks, predictions, intensity);
																								// repeating the same for ActiveLongNew tracks (new tracks added in BackSTB)

			for (deque<Track>::iterator tr = s.bufferTracks.begin(); tr != s.bufferTracks.end(); ) {
				if (!(tr->Exists(prevFrame)) && tr->Exists(currFrame)) {
					int case_num = Predictor(s, prevFrame, tr, unlinkedTracks, predictions, intensity);
					switch (case_num) {
					case 0:
						tr = s.bufferTracks.erase(tr); //No track is connected, it is moved to unlinkedTracks temporarily.
						break;
					case 1:  // the track is connected to a track after this track
						++ tr;
						break;
					default: // the track is connected to a track before this track
						break;
					}
				} else {
					++ tr;
				}
			}

			deque<string> filename;													// connect them using the residual images
																								// loading the actual cameras images at prev frame for connecting links back in time
			for (int i = 0; i < s.ncams; i++) {													// getting the tiff image names
	//			int j = i + 1;
				// exchanging images of camera 3 n 4
	//			if (i == 2)
	//				j = 4;
	//			if (i == 3)
	//				j = 3;
	//			stringstream tiffname; tiffname << s.tiffaddress << "C00" << j << "H001S0001/" << "cam" << j << "frame" << setfill('0') << setw(4) << prevFrame << ".tif";
	//			filename[i] = tiffname.str();
				filename.push_back(s.imgSequence[i][prevFrame - 1]);

				//filename[i] = s.tiffaddress + "frame" + to_string(prevFrame) + "cam" + to_string(i + 1) + ".tif";
			}
			Tiff2DFinder t(s.ncams, s._ipr.Get_threshold(), filename);
			t.FillPixels(pixels_orig);

			Frame tracked;																		// getting the list of tracked particles at prevFrame (from activeLong, activeLongNew and exit tracks)
			for (deque<Track>::iterator tr = s.bufferTracks.begin(); tr != s.bufferTracks.end(); ++tr)
				if (tr->Exists(prevFrame))
					tracked.Add(tr->GetPos(prevFrame - tr->GetTime(0)));


			s._ipr.ReprojImage(tracked, s.OTFcalib, pixels_reproj, STBflag);

			for (int n = 0; n < s.ncams; n++) 													// removing the tracked particles from original image
				for (int i = 0; i < s.Npixh; i++)
					for (int j = 0; j < s.Npixw; j++) {
						int residual = (pixels_orig[n][i][j] - s.fpt*pixels_reproj[n][i][j]);	// using the multiplying factor (fpt) to remove all traces of tracked particles
						pixels_orig[n][i][j] = (residual < 0) ? 0 : residual;
					}
																							// trying to find a link for all the unlinked (exit and activeLong) tracks using residual images
			Shake(s, predictions, intensity, unlinkedTracks);									// shaking the predictions of unlinked tracks

			for (int i = 0; i < predictions.NumParticles(); i++) 								// adding the corrected particle position in prevFrame to its respective track
				unlinkedTracks[i].AddFront(predictions[i], prevFrame);

	//		start2 = clock();
			cout << "Done (" <<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()- start) .count() << "s)" << endl
					<< "\t\tIPR on residuals for new tracks: "<<endl;

			start =std::chrono::system_clock::now();

			if (predictions.NumParticles() != 0) {
				s._ipr.ReprojImage(predictions, s.OTFcalib, pixels_reproj, STBflag);			// updating the original image by removing the newly tracked particles (shaked predictions)
				for (int n = 0; n < s.ncams; n++)
					for (int i = 0; i < s.Npixh; i++)
						for (int j = 0; j < s.Npixw; j++) {
							int residual = (pixels_orig[n][i][j] - s.fpt*pixels_reproj[n][i][j]);
							pixels_orig[n][i][j] = (residual < 0) ? 0 : residual;
						}
			}
																								// applying ipr on the remaining particles in orig images to obtain particle candidates
			s._ipr.DoingSTB(true);
			Frame candidates = s.IPRonResidual(calib, t, pixels_orig, pixels_reproj, pixels_res, predictions);  //Todo: put tracked particles in addition to predictions.
																								// trying to link each activeShortTrack with a particle candidate in prevFrame
			for (deque<Track>::iterator tr = s.activeShortTracks.begin(); tr != s.activeShortTracks.end(); )
				MakeShortLinkResidual_backward(s, prevFrame, candidates, tr, 5);

			cout<<"Time for IPR on residual:"<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()- start) .count() << "s" << endl;

			start =std::chrono::system_clock::now();
			// PRUNING / ARRANGING THE TRACKS
			double thresh = 1.5 * s.largestShift;
			for (deque<Track>::iterator tr = unlinkedTracks.begin(); tr != unlinkedTracks.end(); ) {
				double d1 = pow(Distance(tr->First(), tr->Second()),0.5), d2 = pow(Distance(tr->Second(), tr->Third()),0.5);
				double threshRel = s.maxRelShiftChange*d2, length = tr->Length();
																												// moving all activeLongTracks with displacement more than LargestExp shift to inactiveTracks
				if (d1 > thresh) {
		//					inactiveTracks.push_back(*tr);
					tr->DeleteFront(); //delete the last point
					if (length >= 7) {
						if (tr->Exists(s.last)) {
							s.activeLongTracks.push_back(*tr);
						} else {
							s.inactiveLongTracks.push_back(*tr);
						}
					}
	//				s_al++; a_is++;
				} else if (abs(d1 - d2) > s.maxAbsShiftChange || abs(d1 - d2) > threshRel) { // moving all activeLongTracks with large change in particle shift to inactive/inactiveLong tracks
					tr->DeleteFront();
					if (length >= 7) {
						if (tr->Exists(s.last)) {
							s.activeLongTracks.push_back(*tr);
						} else {
							s.inactiveLongTracks.push_back(*tr);
						}
					}
				} else {
					s.bufferTracks.push_back(*tr); // a new particle is kept and move back to bufferTracks
				}
				tr = unlinkedTracks.erase(tr);
			}
			// moving all activeShortTracks longer than 3 particles to activeLongTracks
			for (deque<Track>::iterator tr = s.activeShortTracks.begin(); tr != s.activeShortTracks.end(); ) {
				if (tr->Length() > 3) {
					s.bufferTracks.push_back(*tr);
					tr = s.activeShortTracks.erase(tr);
				}
					else
					++tr;
			}
			for (Frame::const_iterator pID = candidates.begin(); pID != candidates.end(); ++pID) {	// adding all the untracked candidates to a new short track
				Track startTrack(*pID, prevFrame);
				s.activeShortTracks.push_back(startTrack);
			}

			cout << "Time for pruning:"<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()- start) .count() << "s" << endl;
	//		cout << "Done (" << (clock() - start2) / 1000 << "s)" << endl;

			cout << "\t\tNo. of active Short tracks:	" << s.activeShortTracks.size() << endl;
			cout << "\t\tNo. of active Long tracks:	" << s.activeLongTracks.size() << endl;
			cout << "\t\tNo. of inactive Long tracks:	" << s.inactiveLongTracks.size() << endl;
			cout << "\t\tNo. of buffer tracks:\t\t" << s.bufferTracks.size() << endl;
			cout << "\t\tNo. of exited tracks:		" << s.exitTracks.size() << endl;
	//		cout << "\t\tNo. of inactive tracks:		" << s.inactiveTracks.size() << endl;
	//		cout << "\t\tTime taken for BackSTB at frame " << prevFrame << ": " << (clock() - start0) / 1000 << "s" << endl;


			if (to_save_data) {
				s.MatTracksSave(address, to_string(prevFrame), 1);
			} else {
				//time_t t = time(0);
				if (prevFrame % 100 == 0 || prevFrame == s.first) {  // to debug, every frame should be saved
					cout << "\tSaving the tracks" << endl;

				s.MatTracksSave(address, to_string(prevFrame), 1);
				if (prevFrame % 100 == 0) {
					// remove previous files except anyone that is multiple of 500 for active long track and exit track
					std::remove((address + "ActiveLongTracks" + to_string(prevFrame + 100) + ".txt").c_str());
					std::remove( (address + "ActiveShortTracks" + to_string(prevFrame + 100) + ".txt").c_str());
					std::remove( (address + "BufferTracks" + to_string(prevFrame + 100) + ".txt").c_str());
					// remove previous files except anyone that is multiple of 500 for active long track and exit track
					if ((prevFrame + 100) % 500 != 0) {
						std::remove( (address + "InactiveLongTracks" + to_string(prevFrame + 100) + ".txt").c_str());
						std::remove( (address + "ExitTracks" + to_string(prevFrame + 100) + ".txt").c_str());
					}
				} else if (prevFrame == s.first && prevFrame % 100 != 0) {
					std::remove((address + "ActiveLongTracks" + to_string((int)(prevFrame / 100) * 100 + 100) + ".txt").c_str());
					std::remove( (address + "ActiveShortTracks" + to_string( (int)(prevFrame / 100) * 100 + 100) + ".txt").c_str());
					std::remove( (address + "BufferTracks" + to_string((int)(prevFrame / 100) * 100 + 100) + ".txt").c_str());
					if (((int)(prevFrame / 100) * 100 + 100) % 500 != 0) {
						std::remove( (address + "InactiveLongTracks" + to_string((int)(prevFrame / 100) * 100 + 100) + ".txt").c_str());
						std::remove( (address + "ExitTracks" + to_string((int)(prevFrame / 100) * 100 + 100) + ".txt").c_str());
					}
				}
				if (prevFrame % 500 == 0) {
					// save the inactive long tracks
					string str = to_string(prevFrame);
					string X5 = "InactiveLongTracks" + str;
					s.SaveTrackToTXT(s.inactiveLongTracks, address + X5);
					// empty inactiveLongTracks
					s.inactiveLongTracks.erase(s.inactiveLongTracks.begin(), s.inactiveLongTracks.end());
					// save the inactive long tracks
					string X6 = "ExitTracks" + str;
					s.SaveTrackToTXT(s.exitTracks, address + X6);
					// empty inactiveLongTracks
					s.exitTracks.erase(s.exitTracks.begin(), s.exitTracks.end());
				}

				}
			}
		}

	// moving all bufferTracks to activeLongTracks or inactive long tracks
	for (deque<Track>::iterator tr = s.bufferTracks.begin(); tr != s.bufferTracks.end(); ) {
		if (tr->Exists(s.last)) {
			s.activeLongTracks.push_back(*tr);
		} else {
			s.inactiveLongTracks.push_back(*tr);
		}
		tr = s.bufferTracks.erase(tr);
	}

	s.MatTracksSave(address, to_string(first), 0);

	cout << "\t\tNo. of active Short tracks:	" << s.activeShortTracks.size() << endl;
	cout << "\t\tNo. of active Long tracks:	" << s.activeLongTracks.size() << endl;
	cout << "\t\tNo. of exited tracks:		" << s.exitTracks.size() << endl;
	cout << "\t\tNo. of inactive long tracks:		" << s.inactiveLongTracks.size() << endl;
//	cout << "\tTotal time taken for BackSTB: " << (clock() - start) / 1000 << "s" << endl;
}

int BackSTB::Predictor(STB& s, int prevFrame, deque<Track>::iterator& Ltr,
						deque<Track>& unlinkedTracks, Frame& predictions,
						deque<double>& intensity) {
																								
	double mindistsqr = 1e6;
	deque<deque<Track>::iterator> matches;

	vector< vector<double> > predCoeff(3);													// using polynomial fit / Wiener filter to predict the particle position at prevFrame
	vector<string> direction = { "X", "Y", "Z" };
//	direction[0] = "X"; direction[1] = "Y"; direction[2] = "Z";
	double est[3];
//	for (int i = 0; i < 3; i++) {
//		predCoeff[i] = Polyfit(*Ltr, direction[i], prevFrame);								// predictor coefficients for X->0, Y->1 and Z->2
//		est[i] = predCoeff[i][0] + predCoeff[i][1] * prevFrame + predCoeff[i][2] * pow(prevFrame, 2) + predCoeff[i][3] * pow(prevFrame, 3);
//	}
	for (int i = 0; i < 3; i++) {
		if (Ltr->Length() < 4) {																// if length is 4 or 5, use all points to get 2nd degree polynomial
	//				predCoeff[i] = Polyfit(*tr, direction[i], tr->Length(), 2);
	//				est[i] = predCoeff[i][0] + predCoeff[i][1] * t + predCoeff[i][2] * pow(t, 2);
			est[i] = LMSWienerPredictor(*Ltr, direction[i], 3);
		}
		else if (Ltr->Length() < 6) {														// if length is 6 to 10, use all points to get 3rd degree polynomial
	//				predCoeff[i] = Polyfit(*tr, direction[i], tr->Length(), 3);
	//				//cout << "break" << tr->Length() << endl;
	//				est[i] = predCoeff[i][0] + predCoeff[i][1] * t + predCoeff[i][2] * pow(t, 2) + predCoeff[i][3] * pow(t, 3);
			est[i] = LMSWienerPredictor(*Ltr, direction[i], Ltr->Length() - 1);
		}
		else {																				// if length is more than 11, use last 10 points to get 3rd degree polynomial
	//				predCoeff[i] = Polyfit(*tr, direction[i], 10, 3);
	//				est[i] = predCoeff[i][0] + predCoeff[i][1] * t + predCoeff[i][2] * pow(t, 2) + predCoeff[i][3] * pow(t, 3);
			est[i] = LMSWienerPredictor(*Ltr, direction[i], 5);
		}
	}
	Position estimate(est[0], est[1], est[2]);												// estimated position at prevFrame
																							// checking if it can link back to any inactive track
	for (deque<Track>::iterator Itr = s.bufferTracks.begin(); Itr != s.bufferTracks.end(); ++Itr) {
		if (Itr->Exists(prevFrame)) {
			Position potentialMatch(Itr->GetPos(prevFrame - Itr->GetTime(0)));
			double distsqr = Distance(estimate, potentialMatch);

			if (distsqr < distThreshold) {
				mindistsqr = min(mindistsqr, distsqr);
				if (distsqr <= mindistsqr)
					matches.push_back(Itr);  // if there are several tracks that meet the criteria,
											// Add those which are smaller than the previous tracks to the match
			}
		}
	}

	if (matches.size() > 0) {																// if a link is found in inactiveTracks, 
		int length = matches.back()->Length();
		for (int t = length - 1; t >= 0; t--) {													// delete all points after prevFrame on that inactive track,
			if (matches.back()->GetTime(t) > prevFrame)
				matches.back()->DeleteBack();
			else
				break;
		}

		Ltr->AddFront(*matches.back());														// link the updated inactive track to the beginning of activeTrack and
		s.bufferTracks.erase(matches.back());												// delete the linked track from active long new Tracks
		if (matches.back() > Ltr)
			return 1; // return 1 to indicate the track is linked to another track adter this track
		else
			return 2; // the track is connected to a track before it.
	}
	else {
		// if no link is found,
		if (boundary_check.Check(estimate)) { // if the estimate particle is inside the boundary
			unlinkedTracks.push_back(*Ltr);														// add the unlinked active track and its prediction in prevFrame to a list
			predictions.Add(estimate);															// this list will later be used to find links from the residual images
			intensity.push_back(1);
		} else { //the track should be put into exit tracks since it is outside the boundary
			s.exitTracks.push_back(*Ltr);
		}
		return 0; // return false to indicate no track is connected
	}
}

// predict the next point with Wiener Predictor using LMS algorithm
double BackSTB::LMSWienerPredictor(Track tracks, string direction, int order) {
	int size = tracks.Length();
	double* series = new double[order + 1];

	// get the series to be predicted
	for (unsigned int i = 0; i < order + 1; ++i) {
		if (direction == "X" || direction == "x")
			series[order - i] = tracks[i].X();	//x-values or
		else if (direction == "Y" || direction == "y")
			series[order - i] = tracks[i].Y();	//y-values or
		else if (direction == "Z" || direction == "z")
			series[order - i] = tracks[i].Z();	//z-values
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

// filter for prediction 
// polynomial fit
vector<double> BackSTB::Polyfit(Track tracks, string direction, int prevFrame) {

	int N = 4, size = tracks.Length(), n = 3;
	cout.precision(4);                        //set precision
	cout.setf(ios::fixed);

	double t[4], x[4];

	for (int i = 0; i < N; i++) {
		int startIndex = prevFrame + 1 - tracks.GetTime(0);
		t[N - i - 1] = tracks.GetTime(startIndex + i);	//time

		if (direction == "X" || direction == "x")
			x[N - i - 1] = tracks[startIndex + i].X();	//x-values or
		else if (direction == "Y" || direction == "y")
			x[N - i - 1] = tracks[startIndex + i].Y();	//y-values or
		else if (direction == "Z" || direction == "z")
			x[N - i - 1] = tracks[startIndex + i].Z();	//z-values
	}

	// n = 3 is the degree;
	double T[2 * 3 + 1];                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	for (int i = 0; i < 2 * n + 1; i++) {
		T[i] = 0;
		for (int j = 0; j<N; j++)
			T[i] = T[i] + pow(t[j], i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	}

	double B[3 + 1][3 + 2] = {}, a[3 + 1] = {};            //B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j <= 3; j++)
			B[i][j] = T[i + j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
	}

	double X[3 + 1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	for (int i = 0; i < n + 1; i++) {
		X[i] = 0;
		for (int j = 0; j<N; j++)
			X[i] = X[i] + pow(t[j], i)*x[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	}

	for (int i = 0; i < n + 1; i++)
		B[i][n + 1] = X[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)

	n = n + 1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations
							  //cout << "\nThe Normal(Augmented Matrix) is as follows:\n";

	for (int i = 0; i < n; i++) {                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
		for (int k = i + 1; k < n; k++) {
			if (B[i][i] < B[k][i]) {
				for (int j = 0; j <= n; j++) {
					double temp = B[i][j];
					B[i][j] = B[k][j];
					B[k][j] = temp;
				}
			}
		}
	}

	for (int i = 0; i < n - 1; i++) {           //loop to perform the gauss elimination
		for (int k = i + 1; k < n; k++) {
			double t = B[k][i] / B[i][i];
			for (int j = 0; j <= n; j++) {
				B[k][j] = B[k][j] - t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
			}
		}
	}

	for (int i = n - 1; i >= 0; i--)                //back-substitution
	{                        //x is an array whose values correspond to the values of x,y,z..
		a[i] = B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
		for (int j = 0; j < n; j++) {
			if (j != i)            //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
				a[i] = a[i] - B[i][j] * a[j];
		}
		a[i] = a[i] / B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
	}

	vector<double> coeff;
	for (int i = 0; i < 4; i++) {
		coeff.push_back(a[i]);
	}
	return coeff;
}

// For shaking the predicted estimates
void BackSTB::Shake(STB& s, Frame& estimate, deque<double>& intensity, deque<Track>& unlinkedTracks) {

	if (estimate.NumParticles() > 0) {
		for (int index = 0; index < estimate.NumParticles(); index++) 						// adding 2D image centers and intensity data to the estimates 
			s._ipr.FullData(estimate[index], intensity[index], s.ncams, ALL_CAMS);	
																							// removing ambiguous, ghost and particles that disappear on at least 2 cams
		deque<int> ignoreCam = Rem(s, estimate, intensity, s._ipr.Get_mindist3D(), unlinkedTracks);		

		for (int loopInner = 0; loopInner < s.it_innerloop; loopInner++) {
			double del;
			if (loopInner < 1)  del = config.shaking_shift;
			else if (loopInner < 5)  del = config.shaking_shift / pow(2,loopInner - 1);
			else  del = config.shaking_shift/100;

			s._ipr.ReprojImage(estimate, s.OTFcalib, pixels_reproj, STBflag);				// adding the estimated particles to the reprojected image

			for (int n = 0; n < s.ncams; n++) 												// updating the residual image by removing the estimates
				for (int i = 0; i < s.Npixh; i++)
					for (int j = 0; j < s.Npixw; j++)
						pixels_res[n][i][j] = (pixels_orig[n][i][j] - abs(pixels_reproj[n][i][j]));

#pragma omp parallel num_threads(24)
						{
//			int index = 0;
			// correcting the estimated positions and their intensity by shaking
//			for (Frame::const_iterator pID = estimate.begin(); pID != estimate.end(); ++pID) {
#pragma omp for
			for (int i = 0; i < estimate.NumParticles(); ++i) {
				Frame::const_iterator pID = estimate.begin() + i;
				OTF otf_calib(s.OTFcalib);
				Shaking shake(s.ncams, ignoreCam[i], otf_calib, s.Npixw, s.Npixh, s._ipr.Get_psize(), del, *pID, s.cams, pixels_res, intensity[i]);
				estimate[i] = shake.Get_posnew();
				intensity[i] = shake.Get_int();
//				index++;
			}
						}
		}

		for (int index = 0; index < estimate.NumParticles(); index++) 						// adding 2D image centers and intensity data after shaking
			s._ipr.FullData(estimate[index], intensity[index], s.ncams, ALL_CAMS);
																							
		Rem(s, estimate, intensity, s._ipr.Get_mindist3D(), unlinkedTracks);				// removing ambiguous particles and particles that did not find a match on the actual images	

		s._ipr.ReprojImage(estimate, s.OTFcalib, pixels_reproj, STBflag);

		for (int n = 0; n < ncams; n++) 													// updating the residual image
			for (int i = 0; i < Npixh; i++)
				for (int j = 0; j < Npixw; j++) {
					int residual = (pixels_orig[n][i][j] - s.fpt*pixels_reproj[n][i][j]);	// using the multiplying factor (fpt) to remove all traces of tracked particles
					pixels_orig[n][i][j] = (residual < 0) ? 0 : residual;
				}
	}
}

// removing ghost,ambiguous and particles leaving measurement domain
deque<int> BackSTB::Rem(STB& s, Frame& pos3D, deque<double>& int3D, double mindist_3D, deque<Track>& unlinkedTracks) {

	for (int i = 0; i < pos3D.NumParticles(); i++) 										// deleting particles that are very close to each other
		for (int j = i + 1; j < pos3D.NumParticles(); ) {
			if (Distance(pos3D[i], pos3D[j]) < 4 * mindist_3D * mindist_3D) {
				pos3D.Delete(j); int3D.erase(int3D.begin() + j);
																						// checking if the particle belongs to activeLongNewTrack (track added in BackSTB)
				//if (std::find(s.bufferTracks.begin(), s.bufferTracks.end(), *unlinkedTracks[j]) != s.bufferTracks.end()) {
				//	s.inactiveTracks.push_back(*unlinkedTracks[j]);						// shifting the corresponding activeLongNewTrack to inactiveTracks
				//	s.bufferTracks.erase(unlinkedTracks[j]);
				//}
//				unlinkedTracks.at(j).DeleteFront();
				if (unlinkedTracks.at(j).Length() > 7) {
					if (unlinkedTracks.at(j).Exists(s.last)) {
						s.activeLongTracks.push_back(unlinkedTracks.at(j));
					} else {
						s.inactiveLongTracks.push_back(unlinkedTracks.at(j));
					}
				}
				unlinkedTracks.erase(unlinkedTracks.begin() + j);
			}
			else
				j++;
		}

	int ghost = 0;
//	double avgInt = 0;																	// deleting based on intensity
//	for (int i = 0; i < int3D.size(); i++)
//		if (int3D[i] > 0)
//			avgInt = avgInt + int3D[i];
//
//	avgInt = avgInt / int3D.size();

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

	for (int index = 0; index < int3D.size(); ) {							// removing a particle if its intensity falls below a % of the avg intensity
		if (int3D[index] < s._ipr.Get_intensityLower() * avgInt) {
			pos3D.Delete(index); int3D.erase(int3D.begin() + index);

			// NO NEED TO CHECK, 01.05.2019 Since all unlinked tracks are moved out from the active long new tracks
			// checking if the particle belongs to activeLongNewTrack (track added in BackSTB)
//			//TODO: to check this inorder to remove -fpermissive. Shiyong Tan, 2/1/18
//			deque<Track>::iterator it = *(std::find(&s.bufferTracks.begin(), &s.bufferTracks.end(), unlinkedTracks[index]));
//			if (it != s.bufferTracks.end()) {
//			//	s.inactiveTracks.push_back(*unlinkedTracks[index]);						// shifting the corresponding activeLongNewTrack to inactiveTracks
//			//	s.bufferTracks.erase(unlinkedTracks[index]);
//			}

//			unlinkedTracks.at(index).DeleteFront();
			if (unlinkedTracks.at(index).Length() > 7) {
				if (unlinkedTracks.at(index).Exists(s.last)) {
					s.activeLongTracks.push_back(unlinkedTracks.at(index));
				} else {
					s.inactiveLongTracks.push_back(unlinkedTracks.at(index));
				}
			}

			unlinkedTracks.erase(unlinkedTracks.begin() + index);
			ghost++;
		} else {
			++index;
		}
	}

	deque<int> ignoreCam(pos3D.NumParticles());
	double intensityThresh = 0.25*s._ipr.Get_threshold();

	for (int i = 0; i < pos3D.NumParticles(); ) {										// deleting particles that are outside the image bounds on atleast 2 cameras
		double xlim = (s.Npixh - 1) / 2, ylim = (s.Npixw - 1) / 2;
		int leftcams = 0, belowIntensity = 0; double shiftThreshold = pow(s.largestShift, 2);																	// if the particle disappears only on 1 cam, ignore that cam and perform shaking
		ignoreCam[i] = 100;
		if (abs(pos3D[i].X1() - xlim) > xlim || abs(pos3D[i].Y1() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 0;
		}else if (pixels_orig[0][(int)round(pos3D[i].Y1())][(int)round(pos3D[i].X1())] < intensityThresh) {
					belowIntensity++; ignoreCam[i] = 0;
		}

		if (abs(pos3D[i].X2() - xlim) > xlim || abs(pos3D[i].Y2() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 1;
		}else if (pixels_orig[1][(int)round(pos3D[i].Y2())][(int)round(pos3D[i].X2())] < intensityThresh) {
			belowIntensity++; ignoreCam[i] = 1;
		}

		if (abs(pos3D[i].X3() - xlim) > xlim || abs(pos3D[i].Y3() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 2;
		}else if (pixels_orig[2][(int)round(pos3D[i].Y3())][(int)round(pos3D[i].X3())] < intensityThresh) {
			belowIntensity++; ignoreCam[i] = 2;
		}

		if (abs(pos3D[i].X4() - xlim) > xlim || abs(pos3D[i].Y4() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 3;
		}else if (pixels_orig[3][(int)round(pos3D[i].Y4())][(int)round(pos3D[i].X4())] < intensityThresh) {
			belowIntensity++; ignoreCam[i] = 3;
		}

		if (leftcams >= 2 || belowIntensity >= 2) {															// if the particle disappears on at least 2 cams
			pos3D.Delete(i); int3D.erase(int3D.begin() + i);							// delete the particle 
																						// checking if the particle belongs to activeLongNewTrack (track added in BackSTB)
			//if (std::find(s.bufferTracks.begin(), s.bufferTracks.end(), *unlinkedTracks[i]) != s.bufferTracks.end()) {
			//	s.exitTracks.push_back(*unlinkedTracks[i]);								// shifting the corresponding activeLongNewTrack to exitTracks
			//	s.bufferTracks.erase(unlinkedTracks[i]);
			//}																			// checking if the particle belongs to activeLongTrack 
			//else if (std::find(s.activeLongTracks.begin(), s.activeLongTracks.end(), *unlinkedTracks[i]) != s.activeLongTracks.end()) {
			//	s.exitTracks.push_back(*unlinkedTracks[i]);								// shifting the corresponding activeLongTrack to exitTracks
			//	s.activeLongTracks.erase(unlinkedTracks[i]);
			//}
			if (unlinkedTracks.at(i).Length() > 7) {
				if (unlinkedTracks.at(i).Exists(s.last)) {
					s.activeLongTracks.push_back(unlinkedTracks.at(i));
				} else {
					s.inactiveLongTracks.push_back(unlinkedTracks.at(i));
				}
			}
			unlinkedTracks.erase(unlinkedTracks.begin() + i);
			ignoreCam[i] = 100;
		}
		else
			i++;
	}
	return ignoreCam;
}

// linking short tracks with particle candidates in residual images
void BackSTB::MakeShortLinkResidual_backward(STB& s,int prevFrame, Frame& candidates, deque<Track>::iterator& tr, int iterations) {

	pair<Frame::const_iterator, float> cost;

	for (int it = 0; it < iterations; it++) {													// iteratively trying to find a link for the short track from particle candidates
		double rsqr = pow(pow(1.1, it) * 3 * s.avgIPDist, 2);
		double shift = pow(1.1, it) * s.largestShift;
		deque<double> dist;
		double totalDist = 0;
		deque<Position> disp;
		Position vel(0, 0, 0);// calculating the predictive vel. field as an avg. of particle velocities from neighbouring tracks
		double d;

		for (unsigned int j = 0; j < s.activeLongTracks.size(); j++) {									// identifying the neighbouring tracks (using 3*avg interparticle dist.) and getting their particle velocities
			int currFrameIndex = prevFrame + 1 - s.activeLongTracks[j].GetTime(0);
			if (s.activeLongTracks[j].Exists(prevFrame + 1) && s.activeLongTracks[j].Exists(prevFrame)) {
				double dsqr = Distance(tr->GetPos(0), s.activeLongTracks[j][currFrameIndex]);
				if (dsqr < rsqr) {
					d = pow(dsqr, .5);
					totalDist = totalDist + d;
					dist.push_back(d);
					disp.push_back(s.activeLongTracks[j][currFrameIndex - 1] - s.activeLongTracks[j][currFrameIndex]);
				}
			}		
		}

		for (unsigned int j = 0; j < s.exitTracks.size(); j++) {											// identifying the neighbouring tracks from exit tracks
			int currFrameIndex = prevFrame + 1 - s.exitTracks[j].GetTime(0);
			if (s.exitTracks[j].Exists(prevFrame + 1) && s.exitTracks[j].Exists(prevFrame)) {
				double dsqr = Distance(tr->GetPos(0), s.exitTracks[j][currFrameIndex]);
				if (dsqr < rsqr) {
					d = pow(dsqr, .5);
					totalDist = totalDist + d;
					dist.push_back(d);
					disp.push_back(s.exitTracks[j][currFrameIndex - 1] - s.exitTracks[j][currFrameIndex]);
				}
			}
		}

		for (unsigned int j = 0; j < s.bufferTracks.size(); j++) {											// identifying the neighbouring tracks from exit tracks
			int currFrameIndex = prevFrame + 1 - s.bufferTracks[j].GetTime(0);
			if (s.bufferTracks[j].Exists(prevFrame + 1) && s.bufferTracks[j].Exists(prevFrame)) {
				double dsqr = Distance(tr->GetPos(0), s.bufferTracks[j][currFrameIndex]);
				if (dsqr < rsqr) {
					d = pow(dsqr, .5);
					totalDist = totalDist + d;
					dist.push_back(d);
					disp.push_back(s.bufferTracks[j][currFrameIndex - 1] - s.bufferTracks[j][currFrameIndex]);
				}
			}
		}


		if (dist.size() > 0) {																// if it finds neighbouring tracks
			for (int j = 0; j < dist.size(); j++) {												// perform Gaussian wt. avg. to get velocity field (backward in time)
				double weight = (dist[j] / totalDist);
				vel.Set_X(vel.X() + weight*disp[j].X());
				vel.Set_Y(vel.Y() + weight*disp[j].Y());
				vel.Set_Z(vel.Z() + weight*disp[j].Z());
			}
			Position estimate = tr->Last() + vel;
																								// finding a link for this short track using the velocity field
			cost = s.ComputeCost(candidates, candidates, s.searchRadiusSTB, estimate, vel, vel, true);
		}

		else {																				// if no neighbouring tracks are identified
			Position estimate = tr->GetPos(0);													// using nearest neighbour to make a link
			cost = s.ComputeCost(candidates, candidates, shift, estimate, vel, vel, true);
		}

		if (cost.second != s.UNLINKED) {													// if linked with a candidate
			tr->AddFront(*cost.first, prevFrame);												// add the candidate to the track, delete it from list of untracked candidates, 
			++tr;																				// and break the loop
			candidates.Delete(cost.first.where());
			break;
		}
	}

	if (cost.second == s.UNLINKED) {														// if no link is found for the short track
//		if (tr->Length() > 1) 																	// and if the track has at least 2 particles
//			s.inactiveTracks.push_back(*tr);													// add it to inactive tracks

		tr = s.activeShortTracks.erase(tr);														// then delete from activeShortTracks
	}
}

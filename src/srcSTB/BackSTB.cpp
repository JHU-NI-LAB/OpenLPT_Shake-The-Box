#include <STB.h>
#include <BackSTB.h>

#define ALL_CAMS 100
#define STBflag true
#define IPRflag false

using namespace std;

BackSTB::BackSTB(int firstFrame, int lastFrame, double threshold) : first(firstFrame), last(lastFrame), distThreshold(pow(threshold,2)) {}

void BackSTB::UpdateTracks(STB& s) {
	cout << "\tPerforming BackSTB" << endl;
	s.activeShortTracks.clear();
	clock_t start;
	start = clock();

	ncams = s.ncams;
	Npixh = s.Npixh;
	Npixw = s.Npixw;

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

	for (int currFrame = first + last; currFrame != first; currFrame--) {						// Linking the tracks back in time

		clock_t start0, start1, start2;
		start0 = clock();

		int prevFrame = currFrame - 1;
		cout << "\tBackSTB at frame: " << prevFrame << endl;
		cout << "\t\tCheck Points: " << endl;

		start1 = clock();
		cout << "\t\t\tLinking Long/Exit tracks: ";
		deque<deque<Track>::iterator> unlinkedTracks;
		Frame predictions; deque<double> intensity;						
																							// for each track in activeLongTracks
		for (deque<Track>::iterator Ltr = s.activeLongTracks.begin(); Ltr != s.activeLongTracks.end(); ++Ltr)
			if (!(Ltr->Exists(prevFrame)) && Ltr->Exists(currFrame))					// if the track has no point in the prevFrame and has a point in the currFrame (unlinked track)
				Predictor(s, prevFrame, Ltr, unlinkedTracks, predictions, intensity);		// predict its position in prevFrame (using 4 time steps after prevFrame)

																							// repeating the same for exit tracks
		for (deque<Track>::iterator Etr = s.exitTracks.begin(); Etr != s.exitTracks.end(); ++Etr)
			if (!(Etr->Exists(prevFrame)) && Etr->Exists(currFrame))
				Predictor(s, prevFrame, Etr, unlinkedTracks, predictions, intensity);
																							// repeating the same for ActiveLongNew tracks (new tracks added in BackSTB)
		for (deque<Track>::iterator Ltr = s.activeLongNewTracks.begin(); Ltr != s.activeLongNewTracks.end(); ++Ltr)
			if (!(Ltr->Exists(prevFrame)) && Ltr->Exists(currFrame))					
				Predictor(s, prevFrame, Ltr, unlinkedTracks, predictions, intensity);
		

		deque<string> filename(s.ncams);													// connect them using the residual images
																							// loading the actual cameras images at prev frame for connecting links back in time 
		for (int i = 0; i < s.ncams; i++) {													// getting the tiff image names
			int j = i + 1;
			// exchanging images of camera 3 n 4
			if (i == 2)
				j = 4;
			if (i == 3)
				j = 3;
			stringstream tiffname; tiffname << s.tiffaddress << "C00" << j << "H001S0001/" << "cam" << j << "frame" << setfill('0') << setw(4) << prevFrame << ".tif";
			filename[i] = tiffname.str();
			//filename[i] = s.tiffaddress + "frame" + to_string(prevFrame) + "cam" + to_string(i + 1) + ".tif";
		}
		Tiff2DFinder t(s.ncams, s._ipr.Get_threshold(), filename);
		t.FillPixels(pixels_orig);

		Frame tracked;																		// getting the list of tracked particles at prevFrame (from activeLong, activeLongNew and exit tracks)
		for (deque<Track>::iterator tr = s.activeLongTracks.begin(); tr != s.activeLongTracks.end(); ++tr)
			if (tr->Exists(prevFrame))
				tracked.Add(tr->GetPos(prevFrame - tr->GetTime(0)));

		for (deque<Track>::iterator tr = s.exitTracks.begin(); tr != s.exitTracks.end(); ++tr)
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
			unlinkedTracks[i]->AddFront(predictions[i], prevFrame);

		start2 = clock();
		cout << "Done (" << (clock() - start1) / 1000 << "s)" << endl << "\t\t\tIPR on residuals for new tracks: ";

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
		Frame candidates = s.IPRonResidual(calib, t, pixels_orig, pixels_reproj, pixels_res, predictions);
																							// trying to link each activeShortTrack with a particle candidate in prevFrame
		for (deque<Track>::iterator tr = s.activeShortTracks.begin(); tr != s.activeShortTracks.end(); )
			MakeShortLinkResidual_backward(s, prevFrame, candidates, tr, 5);
																							// moving all activeShortTracks longer than 3 particles to activeLongTracks
		for (deque<Track>::iterator tr = s.activeShortTracks.begin(); tr != s.activeShortTracks.end(); ) {
			if (tr->Length() > 3) {
				s.activeLongNewTracks.push_back(*tr);
				tr = s.activeShortTracks.erase(tr);
			}
			else
				++tr;
		}

		for (Frame::const_iterator pID = candidates.begin(); pID != candidates.end(); ++pID) {	// adding all the untracked candidates to a new short track 	
			Track startTrack(*pID, prevFrame);
			s.activeShortTracks.push_back(startTrack);
		}

		cout << "Done (" << (clock() - start2) / 1000 << "s)" << endl;

		cout << "\t\tNo. of active Short tracks:	" << s.activeShortTracks.size() << endl;
		cout << "\t\tNo. of active Long tracks:	" << s.activeLongTracks.size() << endl;
		cout << "\t\tNo. of active Long new tracks:	" << s.activeLongNewTracks.size() << endl;
		cout << "\t\tNo. of exited tracks:		" << s.exitTracks.size() << endl;
		cout << "\t\tNo. of inactive tracks:		" << s.inactiveTracks.size() << endl;
		cout << "\t\tTime taken for BackSTB at frame " << prevFrame << ": " << (clock() - start0) / 1000 << "s" << endl;
	}

	// moving all activeLongNewTracks to activeLongTracks
	for (deque<Track>::iterator tr = s.activeLongNewTracks.begin(); tr != s.activeLongNewTracks.end(); ) {
		s.activeLongTracks.push_back(*tr);
		tr = s.activeLongNewTracks.erase(tr);
	}

	cout << "\t\tNo. of active Short tracks:	" << s.activeShortTracks.size() << endl;
	cout << "\t\tNo. of active Long tracks:	" << s.activeLongTracks.size() << endl;
	cout << "\t\tNo. of exited tracks:		" << s.exitTracks.size() << endl;
	cout << "\t\tNo. of inactive tracks:		" << s.inactiveTracks.size() << endl;
	cout << "\tTotal time taken for BackSTB: " << (clock() - start) / 1000 << "s" << endl;
}

void BackSTB::Predictor(STB& s, int prevFrame, deque<Track>::iterator& Ltr,
						deque<deque<Track>::iterator>& unlinkedTracks, Frame& predictions, 
						deque<double>& intensity) {
																								
	double mindistsqr = 1e6;
	deque<deque<Track>::iterator> matches;

	vector<vector<double>> predCoeff(3);													// using polynomial fit / Wiener filter to predict the particle position at prevFrame
	vector<string> direction = { "X", "Y", "Z" };
	vector<double> est(3);
	for (int i = 0; i < 3; i++) {
		predCoeff[i] = Polyfit(*Ltr, direction[i], prevFrame);								// predictor coefficients for X->0, Y->1 and Z->2
		est[i] = predCoeff[i][0] + predCoeff[i][1] * prevFrame + predCoeff[i][2] * pow(prevFrame, 2) + predCoeff[i][3] * pow(prevFrame, 3);
	}
	Position estimate(est[0], est[1], est[2]);												// estimated position at prevFrame
/*																							// checking if it can link back to any inactive track
	for (deque<Track>::iterator Itr = s.inactiveTracks.begin(); Itr != s.inactiveTracks.end(); ++Itr) {
		if (Itr->Exists(prevFrame)) {
			Position potentialMatch(Itr->GetPos(prevFrame - Itr->GetTime(0)));
			double distsqr = Distance(estimate, potentialMatch);

			if (distsqr < distThreshold) {
				mindistsqr = min(mindistsqr, distsqr);
				if (distsqr <= mindistsqr)
					matches.push_back(Itr);
			}
		}
	}

	if (matches.size() > 0) {																// if a link is found in inactiveTracks, 
		int length = matches.back()->Length();
		for (int t = length; t >= 0; t--) {													// delete all points after prevFrame on that inactive track,										
			if (matches.back()->GetTime(t) > prevFrame)
				matches.back()->DeleteBack();
			else
				break;
		}

		Ltr->AddFront(*matches.back());														// link the updated inactive track to the beginning of activeTrack and
		s.inactiveTracks.erase(matches.back());												// delete the linked track from inactiveTracks
	}
	else {																					// if no link is found,
*/
		unlinkedTracks.push_back(Ltr);														// add the unlinked active track and its prediction in prevFrame to a list
		predictions.Add(estimate);															// this list will later be used to find links from the residual images	
		intensity.push_back(1);
	//}
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
void BackSTB::Shake(STB& s, Frame& estimate, deque<double>& intensity, deque<deque<Track>::iterator>& unlinkedTracks) {

	if (estimate.NumParticles() > 0) {
		for (int index = 0; index < estimate.NumParticles(); index++) 						// adding 2D image centers and intensity data to the estimates 
			s._ipr.FullData(estimate[index], intensity[index], s.ncams, ALL_CAMS);	
																							// removing ambiguous, ghost and particles that disappear on at least 2 cams
		deque<int> ignoreCam = Rem(s, estimate, intensity, s._ipr.Get_mindist3D(), unlinkedTracks);		

		for (int loopInner = 0; loopInner < s.it_innerloop; loopInner++) {
			double del;
			if (loopInner < 2)  del = 0.01;
			else if (loopInner < 5)  del = 0.001;
			else  del = 0.0001;

			s._ipr.ReprojImage(estimate, s.OTFcalib, pixels_reproj, IPRflag);				// adding the estimated particles to the reprojected image

			for (int n = 0; n < s.ncams; n++) 												// updating the residual image by removing the estimates
				for (int i = 0; i < s.Npixh; i++)
					for (int j = 0; j < s.Npixw; j++)
						pixels_res[n][i][j] = (pixels_orig[n][i][j] - abs(pixels_reproj[n][i][j]));


			int index = 0;																	// correcting the estimated positions and their intensity by shaking
			for (Frame::const_iterator pID = estimate.begin(); pID != estimate.end(); ++pID) {
				Shaking shake(s.ncams, ignoreCam[index], s.OTFcalib, s.Npixw, s.Npixh, s._ipr.Get_psize(), del, *pID, s.cams, pixels_res, intensity[index]);
				estimate[index] = shake.Get_posnew();
				intensity[index] = shake.Get_int();
				index++;
			}
		}

		for (int index = 0; index < estimate.NumParticles(); index++) 						// adding 2D image centers and intensity data after shaking
			s._ipr.FullData(estimate[index], intensity[index], s.ncams, ALL_CAMS);
																							
		Rem(s, estimate, intensity, s._ipr.Get_mindist3D(), unlinkedTracks);				// removing ambiguous particles and particles that did not find a match on the actual images	

		for (int n = 0; n < ncams; n++) 													// updating the residual image
			for (int i = 0; i < Npixh; i++)
				for (int j = 0; j < Npixw; j++) {
					int residual = (pixels_orig[n][i][j] - s.fpt*pixels_reproj[n][i][j]);	// using the multiplying factor (fpt) to remove all traces of tracked particles
					pixels_orig[n][i][j] = (residual < 0) ? 0 : residual;
				}
	}
}

// removing ghost,ambiguous and particles leaving measurement domain
deque<int> BackSTB::Rem(STB& s, Frame& pos3D, deque<double>& int3D, double mindist_3D, deque<deque<Track>::iterator>& unlinkedTracks) {

	for (int i = 0; i < pos3D.NumParticles(); i++) 										// deleting particles that are very close to each other
		for (int j = i + 1; j < pos3D.NumParticles(); ) {
			if (Distance(pos3D[i], pos3D[j]) < mindist_3D * mindist_3D) {
				pos3D.Delete(j); int3D.erase(int3D.begin() + j);
																						// checking if the particle belongs to activeLongNewTrack (track added in BackSTB)
				//if (std::find(s.activeLongNewTracks.begin(), s.activeLongNewTracks.end(), *unlinkedTracks[j]) != s.activeLongNewTracks.end()) {
				//	s.inactiveTracks.push_back(*unlinkedTracks[j]);						// shifting the corresponding activeLongNewTrack to inactiveTracks
				//	s.activeLongNewTracks.erase(unlinkedTracks[j]);
				//}
			unlinkedTracks.erase(unlinkedTracks.begin() + j);
			}
			else
				j++;
		}

	int ghost = 0;
	double avgInt = 0;																	// deleting based on intensity
	for (int i = 0; i < int3D.size(); i++) 
		if (int3D[i] > 0) 
			avgInt = avgInt + int3D[i];

	avgInt = avgInt / int3D.size();

	for (int index = int3D.size() - 1; index >= 0; index--)								// removing a particle if its intensity falls below a % of the avg intensity																				
		if (int3D[index] < s._ipr.Get_intensityLower() * avgInt) {
			pos3D.Delete(index); int3D.erase(int3D.begin() + index);
																						// checking if the particle belongs to activeLongNewTrack (track added in BackSTB)
			//TODO: to check this inorder to remove -fpermissive. Shiyong Tan, 2/1/18
			deque<Track>::iterator it = *(std::find(&s.activeLongNewTracks.begin(), &s.activeLongNewTracks.end(), unlinkedTracks[index]));
			if (it != s.activeLongNewTracks.end()) {
			//	s.inactiveTracks.push_back(*unlinkedTracks[index]);						// shifting the corresponding activeLongNewTrack to inactiveTracks
			//	s.activeLongNewTracks.erase(unlinkedTracks[index]);
			}
			unlinkedTracks.erase(unlinkedTracks.begin() + index);
			ghost++;
		}

	deque<int> ignoreCam(pos3D.NumParticles());

	for (int i = 0; i < pos3D.NumParticles(); ) {										// deleting particles that are outside the image bounds on atleast 2 cameras
		double xlim = (s.Npixh - 1) / 2, ylim = (s.Npixw - 1) / 2;
		int leftcams = 0;																// if the particle disappears only on 1 cam, ignore that cam and perform shaking
		ignoreCam[i] = 100;
		if (abs(pos3D[i].X1() - xlim) > xlim || abs(pos3D[i].Y1() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 0;
		}

		if (abs(pos3D[i].X2() - xlim) > xlim || abs(pos3D[i].Y2() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 1;
		}

		if (abs(pos3D[i].X3() - xlim) > xlim || abs(pos3D[i].Y3() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 2;
		}

		if (abs(pos3D[i].X4() - xlim) > xlim || abs(pos3D[i].Y4() - ylim) > ylim) {
			leftcams++; ignoreCam[i] = 3;
		}

		if (leftcams >= 2) {															// if the particle disappears on at least 2 cams
			pos3D.Delete(i); int3D.erase(int3D.begin() + i);							// delete the particle 
																						// checking if the particle belongs to activeLongNewTrack (track added in BackSTB)
			//if (std::find(s.activeLongNewTracks.begin(), s.activeLongNewTracks.end(), *unlinkedTracks[i]) != s.activeLongNewTracks.end()) {
			//	s.exitTracks.push_back(*unlinkedTracks[i]);								// shifting the corresponding activeLongNewTrack to exitTracks
			//	s.activeLongNewTracks.erase(unlinkedTracks[i]);
			//}																			// checking if the particle belongs to activeLongTrack 
			//else if (std::find(s.activeLongTracks.begin(), s.activeLongTracks.end(), *unlinkedTracks[i]) != s.activeLongTracks.end()) {
			//	s.exitTracks.push_back(*unlinkedTracks[i]);								// shifting the corresponding activeLongTrack to exitTracks
			//	s.activeLongTracks.erase(unlinkedTracks[i]);
			//}

			unlinkedTracks.erase(unlinkedTracks.begin() + i);
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
		Position vel(0, 0, 0);																	// calculating the predictive vel. field as an avg. of particle velocities from neighbouring tracks

		for (int j = 0; j < s.activeLongTracks.size(); j++) {									// identifying the neighbouring tracks (using 3*avg interparticle dist.) and getting their particle velocities
			int currFrameIndex = prevFrame + 1 - s.activeLongTracks[j].GetTime(0);
			if (s.activeLongTracks[j].Exists(prevFrame + 1) && s.activeLongTracks[j].Exists(prevFrame)) {
				double dsqr = Distance(tr->GetPos(0), s.activeLongTracks[j][currFrameIndex]);
				if (dsqr < rsqr) {
					totalDist = totalDist + dsqr;
					dist.push_back(dsqr);
					disp.push_back(s.activeLongTracks[j][currFrameIndex - 1] - s.activeLongTracks[j][currFrameIndex]);
				}
			}		
		}

		for (int j = 0; j < s.exitTracks.size(); j++) {											// identifying the neighbouring tracks from exit tracks
			int currFrameIndex = prevFrame + 1 - s.exitTracks[j].GetTime(0);
			if (s.exitTracks[j].Exists(prevFrame + 1) && s.exitTracks[j].Exists(prevFrame)) {
				double dsqr = Distance(tr->GetPos(0), s.exitTracks[j][currFrameIndex]);
				if (dsqr < rsqr) {
					totalDist = totalDist + dsqr;
					dist.push_back(dsqr);
					disp.push_back(s.exitTracks[j][currFrameIndex - 1] - s.exitTracks[j][currFrameIndex]);
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
		if (tr->Length() > 1) 																	// and if the track has at least 2 particles
			s.inactiveTracks.push_back(*tr);													// add it to inactive tracks

		tr = s.activeShortTracks.erase(tr);														// then delete from activeShortTracks
	}
}

#include <STB.h>
#include <ForwardSTB.h>

#define ALL_CAMS 100
#define STBflag true
#define IPRflag false

using namespace std;

ForwardSTB::ForwardSTB(int firstFrame, int lastFrame, double threshold) : first(firstFrame), last(lastFrame), distThreshold(pow(threshold, 2)) {}

void ForwardSTB::UpdateTracks(STB& s) {
	s.activeShortTracks.clear();
	cout << "\tPerforming ForwardSTB" << endl;
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
			pixels_orig[n][i] = new int[s.Npixw]{};
			pixels_reproj[n][i] = new int[s.Npixw]{};
			pixels_res[n][i] = new int[s.Npixw]{};
		}
	}

	Calibration calib(s._ipr.Get_calibfile(), ALL_CAMS, s._ipr.Get_mindist2D(), s._ipr.Get_mindist3D(), ncams);

	for (int currFrame = first - 1; currFrame != last; currFrame++) {						// Linking the tracks back in time
			cout << "\tPerforming BackSTB" << endl;
		clock_t start0, start1, start2;
		start0 = clock();

		int nextFrame = currFrame + 1;
		cout << "\tForwardSTB at frame: " << nextFrame << endl;
		cout << "\t\tCheck Points: " << endl;

		start1 = clock();
		cout << "\t\t\tLinking Long/Exit tracks: ";
		deque<deque<Track>::iterator> unlinkedTracks;
		Frame predictions; deque<double> intensity;
		// for each track in activeLongTracks
		for (deque<Track>::iterator Ltr = s.activeLongTracks.begin(); Ltr != s.activeLongTracks.end(); ++Ltr)
			if (!(Ltr->Exists(nextFrame)) && Ltr->Exists(currFrame) && Ltr->Exists(currFrame - 3))	// if the track has no point in the nextFrame and has a point in the currFrame (unlinked track)
				Predictor(s, nextFrame, Ltr, unlinkedTracks, predictions, intensity);				// predict its position in nextFrame (using 4 time steps before nextFrame)

																									// repeating the same for exit tracks
		for (deque<Track>::iterator Etr = s.exitTracks.begin(); Etr != s.exitTracks.end(); ++Etr)
			if (!(Etr->Exists(nextFrame)) && Etr->Exists(currFrame) && Etr->Exists(currFrame - 3))
				Predictor(s, nextFrame, Etr, unlinkedTracks, predictions, intensity);


		deque<string> filename(s.ncams);															// connect them using the residual images
																									// loading the actual cameras images at nextFrame for connecting links forward in time 
		for (int i = 0; i < s.ncams; i++) 															// getting the tiff image names
			filename[i] = s.tiffaddress + "frame" + to_string(nextFrame) + "cam" + to_string(i + 1) + ".tif";

		Tiff2DFinder t(s.ncams, s._ipr.Get_threshold(), filename);
		t.FillPixels(pixels_orig);

		Frame tracked;																				// getting the list of tracked particles at nextFrame (from activeLong and exit tracks)
		for (deque<Track>::iterator tr = s.activeLongTracks.begin(); tr != s.activeLongTracks.end(); ++tr)
			if (tr->Exists(nextFrame))
				tracked.Add(tr->GetPos(nextFrame - tr->GetTime(0)));

		for (deque<Track>::iterator tr = s.exitTracks.begin(); tr != s.exitTracks.end(); ++tr)
			if (tr->Exists(nextFrame))
				tracked.Add(tr->GetPos(nextFrame - tr->GetTime(0)));

		s._ipr.ReprojImage(tracked, s.OTFcalib, pixels_reproj, STBflag);

		for (int n = 0; n < s.ncams; n++) 															// removing the tracked particles from original image
			for (int i = 0; i < s.Npixh; i++)
				for (int j = 0; j < s.Npixw; j++) {
					int residual = (pixels_orig[n][i][j] - s.fpt*pixels_reproj[n][i][j]);			// using the multiplying factor (fpt) to remove all traces of tracked particles
					pixels_orig[n][i][j] = (residual < 0) ? 0 : residual;
				}
																									// trying to find a link for all the unlinked (exit and activeLong) tracks using residual images
		Shake(s, predictions, intensity, unlinkedTracks);											// shaking the predictions of unlinked tracks

		for (int i = 0; i < predictions.NumParticles(); i++) 										// adding the corrected particle position in nextFrame to its respective track 
			unlinkedTracks[i]->AddNext(predictions[i], nextFrame);

		start2 = clock();
		cout << "Done (" << (clock() - start1) / 1000 << "s)" << endl << "\t\t\tIPR on residuals for new tracks: ";

		if (predictions.NumParticles() != 0) {
			s._ipr.ReprojImage(predictions, s.OTFcalib, pixels_reproj, STBflag);					// updating the original image by removing the newly tracked particles (shaked predictions)
			for (int n = 0; n < s.ncams; n++)
				for (int i = 0; i < s.Npixh; i++)
					for (int j = 0; j < s.Npixw; j++) {
						int residual = (pixels_orig[n][i][j] - s.fpt*pixels_reproj[n][i][j]);
						pixels_orig[n][i][j] = (residual < 0) ? 0 : residual;
					}
		}
																									// applying ipr on the remaining particles in orig images to obtain particle candidates
		Frame candidates = s.IPRonResidual(calib, t, pixels_orig, pixels_reproj, pixels_res, predictions);
																									// trying to link each activeShortTrack with a particle candidate in nextFrame
		for (deque<Track>::iterator tr = s.activeShortTracks.begin(); tr != s.activeShortTracks.end(); )
			s.MakeShortLinkResidual(nextFrame, candidates, tr, 5);
																									// moving all activeShortTracks longer than 3 particles to activeLongTracks
		for (deque<Track>::iterator tr = s.activeShortTracks.begin(); tr != s.activeShortTracks.end(); ) {
			if (tr->Length() > 3) {
				s.activeLongTracks.push_back(*tr);
				tr = s.activeShortTracks.erase(tr);
			}
			else
				++tr;
		}

		for (Frame::const_iterator pID = candidates.begin(); pID != candidates.end(); ++pID) {		// adding all the untracked candidates to a new short track 	
			Track startTrack(*pID, nextFrame);
			s.activeShortTracks.push_back(startTrack);
		}

		cout << "Done (" << (clock() - start2) / 1000 << "s)" << endl;

		cout << "\t\tNo. of active Short tracks:	" << s.activeShortTracks.size() << endl;
		cout << "\t\tNo. of active Long tracks:	" << s.activeLongTracks.size() << endl;
		cout << "\t\tNo. of exited tracks:		" << s.exitTracks.size() << endl;
		cout << "\t\tNo. of inactive tracks:		" << s.inactiveTracks.size() << endl;
		cout << "\t\tTime taken for ForwardSTB at frame " << nextFrame << ": " << (clock() - start0) / 1000 << "s" << endl;
	}

	cout << "\tTotal time taken for ForwardSTB: " << (clock() - start) / 1000 << "s" << endl;
}

void ForwardSTB::Predictor(STB& s, int nextFrame, deque<Track>::iterator& Ltr,
	deque<deque<Track>::iterator>& unlinkedTracks, Frame& predictions,
	deque<double>& intensity) {

	vector< vector<double> > predCoeff(3);													// using polynomial fit / Wiener filter to predict the particle position at nextFrame
	vector<string> direction;
	direction[0] = "X"; direction[1] = "Y"; direction[2] = "Z";
	vector<double> est(3);
	for (int i = 0; i < 3; i++) {
		predCoeff[i] = s.Polyfit(*Ltr, direction[i],4,3);										// predictor coefficients for X->0, Y->1 and Z->2
		est[i] = predCoeff[i][0] + predCoeff[i][1] * nextFrame + predCoeff[i][2] * pow(nextFrame, 2) + predCoeff[i][3] * pow(nextFrame, 3);
	}
	Position estimate(est[0], est[1], est[2]);												// estimated position at nextFrame
																							
	unlinkedTracks.push_back(Ltr);															// add the unlinked active track and its prediction in nextFrame to a list
	predictions.Add(estimate);																// this list will later be used to find links from the residual images	
	intensity.push_back(1);
}

// For shaking the predicted estimates
void ForwardSTB::Shake(STB& s, Frame& estimate, deque<double>& intensity, deque<deque<Track>::iterator>& unlinkedTracks) {

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
deque<int> ForwardSTB::Rem(STB& s, Frame& pos3D, deque<double>& int3D, double mindist_3D, deque<deque<Track>::iterator>& unlinkedTracks) {

	for (int i = 0; i < pos3D.NumParticles(); i++) 										// deleting particles that are very close to each other
		for (int j = i + 1; j < pos3D.NumParticles(); ) {
			if (Distance(pos3D[i], pos3D[j]) < mindist_3D * mindist_3D) {
				pos3D.Delete(j); int3D.erase(int3D.begin() + j);
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
			unlinkedTracks.erase(unlinkedTracks.begin() + i);
		}
		else
			i++;
	}
	return ignoreCam;
}

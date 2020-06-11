#include <NumDataIO.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <ratio>
#include <chrono>

#include "IPR.h"
#include "Common.h"

using namespace std;
#define ALL_CAMS 100
#define REDUCED_CAMS 0
#define STBflag true
#define IPRflag false

IPR::IPR(string& fname, int ncams) : ncams(ncams)
{
	/*
	File contains
	--------------------------------
	# Path to cam calib file
	# Path to tiff files
	# Path to OTF .mat file
	# avg. particle size in pixels_orig (in 1D)
	# No. of outerloop iterations
	# No. of innerloop iterations
	# threshold for 2D finder
	# no. of bits/pixel
	# use reduce cams yes / no
	# no. of loops for each reduced camera combination
	-------------------------------
	*/
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

	int tri; triangulationOnly = false;									// Only use triangulation?
	parsed >> tri;
	if (tri != 0) 
		triangulationOnly = true;

	int ipr; IPROnly = false;
	parsed >> ipr;
	if (ipr != 0)
		IPROnly = true;

	parsed >> calibfile;												// taking the path to calib file
	parsed >> tiffaddress;												// taking the path to tiff files
	parsed >> otfFile;													// taking the path to .mat / .txt file with OTF parameters
	parsed >> psize;													// particle size in x or y direction

	parsed >> it_outerloop;												// # of outer and innerloop iterations
	parsed >> it_innerloop;
	cout << "# of outerloops: " << it_outerloop << endl;
	cout << "# of innerloops: " << it_innerloop << endl;

	parsed >> threshold;												// threshold for 2D finder
	int bitMax;															// max pixel value
	parsed >> bitMax;
	pixelMax = (int)pow(2, bitMax) - 1;

	parsed >> intensityLower;											// ghost particle intensity threshold
	parsed >> mindist_2D;												// reading in tolerences
	parsed >> mindist_3D;
	mindist_2D = mindist_2D * config.factor; // unit conversion
	mindist_3D = mindist_3D * config.factor;

	int reduced;														// reduced cams?
	parsed >> reduced;
	if (reduced == 1) 
		reducedCams = true;
	
	parsed >> it_reducedCam;											// reading in tolerences for reduced cams
	parsed >> mindist_2Dr;
	parsed >> mindist_3Dr;
	mindist_2Dr = mindist_2Dr * config.factor;
	mindist_3Dr = mindist_3Dr * config.factor;
	
	ifstream infile2(calibfile.c_str(), ios::in);						// getting camera parameters
	string line2;
	parsed.str("");
	while (getline(infile2, line2)) {
		size_t commentpos = line2.find('#');
		if (commentpos > 0) {
			if (commentpos < string::npos) {
				line2.erase(commentpos);
			}
			parsed << line2 << '\t';
		}
	}
	infile.close();
	int dummy;
	parsed >> dummy;

	for (int i = 0; i < ncams; ++i)										// now read in the parameters for each camera
		camsAll.push_back(Camera(parsed));

	if (ncams >= 4)														// if there are more than 4 cams, reading the first 4 cams into cams4  TODO: Why?
		for (int i = 0; i < 4; ++i)
			cams4.push_back(camsAll[i]);
	else
		for (int i = 0; i < ncams; ++i)
			cams4.push_back(camsAll[i]);

	Npixh = camsAll[0].Get_Npixh();										// # of pixels in each dimension
	Npixw = camsAll[0].Get_Npixw();

	for (int n = 0; n < ncams; n++) {									// initializing original, residual and reprojected matrices
		pixels_orig.push_back(new int*[Npixh]);
		pixels_reproj.push_back(new int*[Npixh]);
		pixels_res.push_back(new int*[Npixh]);
		for (int i = 0; i < Npixh; i++) {
			pixels_orig[n][i] = new int[Npixw] {};
			pixels_reproj[n][i] = new int[Npixw] {};
			pixels_res[n][i] = new int[Npixw] {};
		}
	}
	m_particle_position_addr = tiffaddress + "ParticlePositions/";
	m_doing_STB = false;
	m_reduce_cam_begin = 0;
}

Frame IPR::FindPos3D(deque< deque<string> > imgNames, int frameNumber)  {

	// clearing all the ipr variables
	pos3D.clear();
	ghost3D.clear();
	intensity3D.clear();
	/*
	 * Modified by Shiyong Tan, 3/30/2018
	 * iframes should be cleared before writing new data into it.
	 * Start:
	 */
	iframes.clear();
	// End
	
	frame = frameNumber;
	// creating a filename with all tiff image names at a particular frame
	/* make sure that the images are in the same sequence as the camera numbers*/
	for (int i = 0; i < ncams; i++) 
		filename.push_back(imgNames[i][frame-1]);
	
// converting the images to 2D dynamic array and finding 2D particle centers
	Tiff2DFinder t(ncams, threshold, filename);
	t.FillPixels(pixels_orig);	// pixels_orig will be filled with pixel intensities for all cameras
	for (int camID = 0; camID < ncams; camID++) {		// filling iframes with the 2D positions on each camera
			try {
				ParticleFinder p(pixels_orig[camID], Npixh, Npixw);//, t.Get_colors(), threshold);
				if (debug_mode == SKIP_IPR_2D_POSITION && frame - 1 < debug_frame_number) { // read 2D position directly
					iframes.push_back(p.ReadParticle2DCenter(imgNames[camID][frame - 1]));
					if (error == NO_FILE) {
						cout<<"The file for reading particle 2D center can't be opened!";
						error = 0;
						goto Find2DCenter;
					}
				} else {
					Find2DCenter:
					p.GetParticle2DCenter(t.Get_colors(), threshold);
					iframes.push_back(p.CreateFrame());
					p.SaveParticle2DCenter(imgNames[camID][frame - 1]);
				}
			}
			catch (out_of_range& e) {
				cerr << e.what() << endl;
				throw runtime_error("Caught out_of_range in ImageSequence::Particle2DList()");
			}
	}

	//Load_2Dpoints("S:/Projects/Bubble/09.28.17/Bubbles and Particles - 250fps/frame", frame, ALL_CAMS);

	int ncams4 = (ncams > 4) ? 4 : ncams;  //No idea why it was set like this.
//	int ncams4 = ncams;
	Calibration calib(cams4, mindist_2D, mindist_3D, ncams4);
	OTF OTFcalib(ncams, otfFile);

	double del = 0.01; //mm in 3D

	clock_t start0;
 	start0 = clock();

// iteratively correcting the position of each particle (SHAKING)

// ALL (OR FIRST 4) CAMERA OUTERLOOP STARTS HERE 
	deque<int> camNums;
	for (int i = 0; i < ncams4; i++)
		camNums.push_back(i);

	deque<Camera> cam_param = calib.Get_cam();

	double pix_length = (cam_param[0].Get_wpix() + cam_param[0].Get_hpix()) / 2; // this will be used as the shift increment

	for (int loopOuter = 0; loopOuter < it_outerloop; loopOuter++) {

		// time
		clock_t start;
		start = clock();

		// increasing the 2D threshold by 10% in each iteration
//		calib.Set_min2D(pow(1.1, loopOuter)*mindist_2D);
		calib.Set_min2D(mindist_2D + pix_length / 4 * loopOuter);
		
		// update pos3Dnew with 3D particles from IPR of current outerloop
		Frame pos3Dnew = IPRLoop(calib, OTFcalib, camNums, ALL_CAMS, t.Get_colors(), pixels_orig, pixels_reproj, pixels_res, loopOuter);
		
		// add the current outerloop particles to pos3D
		pos3D.insert(pos3D.end(), pos3Dnew.begin(), pos3Dnew.end());

		double duration = clock() - start;
		cout << "\t# of particles detected in outerloop" << loopOuter << ": " << pos3Dnew.NumParticles() << endl;
		cout << "\t Time taken for outerloop " << loopOuter << ": " << duration/(double) CLOCKS_PER_SEC << "s" << endl;
		cout << "\tTotal particles (" << pos3D.size() << ")" << endl;

	} 
	
// REDUCED CAMERA LOOP (ignoring 1 camera) 
	if (reducedCams) {
		m_reduce_cam_begin = 1;
		cout << "\t\t\t\t Applying IPR for reduced cams" << endl;

		Calibration calibReduced(calibfile, REDUCED_CAMS, mindist_2Dr, mindist_3Dr, ncams4); // new calibration file for reduced cameras	


		for (int loopOuter = 0; loopOuter < it_reducedCam; loopOuter++) {
			// increasing the 2D threshold by 10% in each iteration
//			calibReduced.Set_min2D(pow(1.1,loopOuter)*mindist_2Dr);
			calibReduced.Set_min2D(mindist_2Dr + pix_length / 4 * loopOuter);
	
			// running IPR by ignoring cameras one by one
			for (int ignoreCam = 0; ignoreCam < ncams4; ignoreCam++) {
				Frame pos3Dnew = IPRLoop(calibReduced, OTFcalib, camNums, ignoreCam, t.Get_colors(), pixels_orig, pixels_reproj, pixels_res, loopOuter);
				cout << "\t # of particles detected in outerloop" << loopOuter << ", ignoring cam" << ignoreCam << ": " << pos3Dnew.NumParticles() << endl;

				pos3D.insert(pos3D.end(), pos3Dnew.begin(), pos3Dnew.end());
			}

			cout << "\t Total particles (" << pos3D.size() << ")" << endl;
		}
		m_reduce_cam_begin = 0;
	}

	double duration0 = clock() - start0;
	cout << "\t Total # of particles detected: " << pos3D.size() << endl;
	cout << "\t Total time taken by IPR: " << duration0/(double) CLOCKS_PER_SEC << "s" << endl;

// FOR TESTING 
	// saving pos3D as a .mat file
	if (triangulationOnly || IPROnly) {
		string reproj = "ReprojectedImage", res = "ResidualImage"+ to_string(frame);
		stringstream mat_name1; mat_name1 << tiffaddress << "posframe" << frameNumber;
		stringstream mat_name2; mat_name2 << tiffaddress << "cleanlistpos2Dframe" << frameNumber;
		SaveParticlePositions(pos3D, mat_name1.str());
		SaveParticlePositions(calib.good2Dpos, mat_name2.str());
		//MatfilePositions(ghost3D, "ghost3D");


	// saving the final residual and reprojected images as .mat files
	/*ReprojImage(pos3D, OTFcalib);
	// shifting the pixels of reprojected image: 1,1 --> 0,0
	for (int n = 0; n < ncams; n++) {
		for (int i = 1; i < Npixh; i++) {
			for (int j = 1; j < Npixw; j++) {
				pixels_reproj[n][i - 1][j - 1] = pixels_reproj[n][i][j];
			}
		}
	}*/
	// shifting the pixels of residual image: 1,1 --> 0,0
	for (int n = 0; n < ncams; n++) {
		for (int i = 1; i < Npixh; i++) {
			for (int j = 1; j < Npixw; j++) {
				pixels_res[n][i - 1][j - 1] = pixels_res[n][i][j];
			}
		}
	}

	//MatfileImage(pixels_reproj, reproj);
	//MatfileImage(pixels_res, res);
	}
	return pos3D;
}


//############################################################################# SUB-FUNCTIONS ##########################################################################

// a function that gives the 3D position by applying a single outerloop of IPR
Frame IPR::IPRLoop(Calibration& calib, OTF& OTFcalib,  deque<int> camNums, int ignoreCam, double colors,
					deque<int**>& orig, deque<int**>& reproj, deque<int**>& res, int loop_time) {
	Frame pos3Dnew;
	double del;
	int id = 0;
	// getting 2D particles from the updated original image (residual after removing identified particles) on reduced cams
	// for the first loop we should not clear iframes
	// for the following loop of IPR in initial phase as well as every loop in the convergence phase, 
	// it should be refreshed by using residual images.
	if (loop_time >= 1 || m_reduce_cam_begin || m_doing_STB) {
//		cout<<(m_reduce_cam_begin)<<endl;
	iframes.clear();
	for (int camID = 0; camID < camNums.size(); camID++)
		if (camNums[camID] != ignoreCam) {
			try {
				ParticleFinder p(orig[camNums[camID]], Npixh, Npixw);//, colors, threshold);
				p.GetParticle2DCenter(colors, threshold);
				iframes.push_back(p.CreateFrame());
				//p.SaveParticle2DCenter("/home/tanshiyong/Documents/Data/Single-Phase/11.03.17/Run1/frame100_" + to_string(camID) + ".txt");
			}
			catch (out_of_range& e) {
				cerr << e.what() << endl;
				throw runtime_error("Caught out_of_range in ImageSequence::Particle2DList()");
			}
		}
	}


//	Load_2Dpoints("S:/Projects/Bubble/Cam_Config_of_10.22.17/10.29.17/BubblesNParticlesHigh_4000fps/BubblesNParicleswithBreakup/Bubble_Reconstruction_Corrected/Bubble_2D_centers", frame, ignoreCam);
	cout <<"\tThe 2D particle number in a image is: "<< iframes[0].NumParticles() << endl;
//	cout << iframes[1].NumParticles() << endl;
//	cout << iframes[2].NumParticles() << endl;
//	cout << iframes[3].NumParticles() << endl;
	// stereomatching to give 3D positions from triangulation
	if ((debug_mode == SKIP_IPR_TRIANGULATION || debug_mode == SKIP_IPR_SHAKING)
			&& frame - 1< debug_frame_number) {
		if (!m_reduce_cam_begin) {
			pos3Dnew = ReadParticlePositions(m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + ".txt");
		} else {
			pos3Dnew = ReadParticlePositions(m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam)
					+ "Loop" + to_string(loop_time) +  ".txt");
		}
		if (error == NO_FILE) {
			string message;
			if (!m_reduce_cam_begin) {
				message = "The triangulation file for frame" + to_string(frame) + "Loop" + to_string(loop_time) + "can't be opened!\n";
			} else {
				message = "The triangulation file for frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam)
							+ "Loop" + to_string(loop_time) + "can't be opened!\n";
			}
			cout<<message;
			error = 0;
			goto StereoMatch;
		} else {
			printf("Read in %i triangulated particles\n", pos3Dnew.NumParticles());
		}
	} else {
		StereoMatch:
		pos3Dnew = calib.Stereomatch(iframes, frame, ignoreCam);
//		SaveParticlePositions(pos3Dnew.Get_PosDeque(), m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + ".txt");
		if (frame < 4) {  // we save the data only for the first 4 frames.
			if (!m_reduce_cam_begin ) {
				SaveParticlePositions(pos3Dnew.Get_PosDeque(), m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + ".txt");
			} else {
				SaveParticlePositions(pos3Dnew.Get_PosDeque(),
						m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) +  ".txt");
			}
		}

	}

	if (((!triangulationOnly) && IPROnly) || !(triangulationOnly || IPROnly)) {
		// initializing the 3D intensity
		deque<double> intensity3Dnew;
		if (debug_mode == SKIP_IPR_SHAKING && frame - 1 < debug_frame_number) {
			string file_path, file_path1;
			if (!m_reduce_cam_begin ) {
				file_path =  m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "Shaking.txt";
				file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "ShakingIntensity.txt";
			} else {
				file_path = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) +  "Shaking.txt";
				file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) +  "ShakingIntensity.txt";
			}
			pos3Dnew = ReadParticlePositions(file_path);
			intensity3Dnew = ReadParticleIntensity(file_path1);
			if (error == NO_FILE) {
				string message;
				if (!m_reduce_cam_begin) {
					message = "The shaking file for frame" + to_string(frame) + "Loop" + to_string(loop_time) + "can't be opened!\n";
				} else {
					message = "The shaking file for frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam)
										+ "Loop" + to_string(loop_time) + "can't be opened!\n";
				}
				cout<<message;
				error = 0;
				goto Shaking;
			} else {
				printf("Read in %i Shaked particles\n", pos3Dnew.NumParticles());
			}
		} else {
			Shaking:
			for (unsigned int i = 0; i < pos3Dnew.NumParticles(); i++)
						intensity3Dnew.push_back(1.0);

					// ######## Innerloop starts here ########
					for (int loopInner = 0; loopInner < it_innerloop; loopInner++) {

						//pair<int, int> bad = Rem(pos3Dnew, intensity3Dnew, mindist_3D);
//						double start = clock();

						if (pos3Dnew.NumParticles() == 0)
							break;

						// Creating the reprojected images (pixels_reproj) by reprojecting the 3D particles onto all cameras using Gaussian ellipse.
						ReprojImage(pos3Dnew, OTFcalib, reproj, IPRflag);
//						ReprojImage(pos3Dnew, OTFcalib, reproj, 0.5);

						// residual image
//						NumDataIO<int> data_io;
//						int* res_pixel = new int[Npixh * Npixw];
//						int* orig_pixel = new int[Npixh * Npixw];
//						int* proj_pixel = new int[Npixh * Npixw];
						for (int n = 0; n < camNums.size(); n++) {
							for (int i = 0; i < Npixh; i++) {
								for (int j = 0; j < Npixw; j++) {
									int residual = (orig[camNums[n]][i][j] - reproj[camNums[n]][i][j]);
									res[camNums[n]][i][j] = residual; //(residual < 0) ? 0 : residual;
//									res_pixel[i * Npixw + j] = residual;
//									orig_pixel[i * Npixw + j] = orig[camNums[n]][i][j];
//									proj_pixel[i * Npixw + j] = reproj[camNums[n]][i][j];
								}
							}
//							data_io.SetTotalNumber(Npixh * Npixw);
//							data_io.SetFilePath("/home/tanshiyong/Documents/Data/SinglePhase/SD78500/ParticlePositions/res.txt");
//							data_io.WriteData(res_pixel);
//							data_io.SetFilePath("/home/tanshiyong/Documents/Data/SinglePhase/SD78500/ParticlePositions/orig.txt");
//							data_io.WriteData(orig_pixel);
//							data_io.SetFilePath("/home/tanshiyong/Documents/Data/SinglePhase/SD78500/ParticlePositions/proj.txt");
//							data_io.WriteData(proj_pixel);
						}

//						if (ignoreCam < ncams) { // add the ignorecam of res 10.19.18 TODO
//							for (int i = 0; i < Npixh; i++) {
//								for (int j = 0; j < Npixw; j++) {
//									int residual = (orig[ignoreCam][i][j] - reproj[ignoreCam][i][j]);
//									res[ignoreCam][i][j] = residual;// (residual < 0) ? 0 : residual;
//									}
//							}
//
//						}


//						delete[] res_pixel;
//						delete[] orig_pixel;
//						delete[] proj_pixel;


						// updating the 3D position and intensity field by shaking
						Frame::const_iterator pIDend = pos3Dnew.end();
//						int index = 0;

						// shaking
//						Frame::const_iterator pID = pos3Dnew.begin();
//						auto start = std::chrono::system_clock::now();
#pragma omp parallel //num_threads(8)
						{//2 4 6 8 10 20 15
//							int TID = omp_get_thread_num();
//							printf("Thread %d is runing\n", TID);
#pragma omp for
						for (int i = 0; i < pos3Dnew.NumParticles(); ++i) {
//						for (Frame::const_iterator pID = pos3Dnew.begin(); pID != pIDend; ++pID) {
							Frame::const_iterator pID = pos3Dnew.begin() + i;
							if (loopInner < 2)  del = mindist_2D;
							else if (loopInner < 5)  del = mindist_2D / pow(2,loopInner - 1);
							else  del = mindist_2D / 20;
//							if (loopInner < 1)  del = config.shaking_shift;			// initial shake TODO
//							else if (loopInner < 5)  del = config.shaking_shift / pow(2,loopInner - 1);//_ipr.mindist_2D/10;	// normal shakes TODO
//							else  del = config.shaking_shift/100;

							OTF otf_calib(OTFcalib);

							Shaking s(ncams, ignoreCam, otf_calib, Npixw, Npixh, psize, del, *pID, camsAll, res, intensity3Dnew[i]);

							pos3Dnew[i] = s.Get_posnew();
							intensity3Dnew[i] = s.Get_int();
//							index++;
//							pID++;
						}
						}
//						double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
//						cout<<"Time taken for shaking in each loop:"<<duration<<endl;
//						auto end = std::chrono::system_clock::now();
//						auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//						cout << elapsed.count() << '\n';

					} // ############# Innerloop ends here #############
					// Save the data
					if (to_save_data) {
					string file_path, file_path1;
						if (!m_reduce_cam_begin ) {
							file_path =  m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "Shaking.txt";
							file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + "ShakingIntensity.txt";
						} else {
							file_path = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) +  "Shaking.txt";
							file_path1 = m_particle_position_addr + "frame" + to_string(frame) + "ReduceCam" + to_string(ignoreCam) + "Loop" + to_string(loop_time) +  "ShakingIntensity.txt";
						}
						SaveParticlePositions(pos3Dnew.Get_PosDeque(),file_path);
						SaveParticleIntensity(intensity3Dnew, file_path1);
					}
		}

		if (pos3Dnew.NumParticles() != 0) {
//			cout <<  "Number of particles after shaking:" <<pos3Dnew.NumParticles()<<endl;
			// save the 3D position, 3D intensity and their correspoding 2D positions
			for (int index = 0; index < pos3Dnew.NumParticles(); index++) 
				FullData(pos3Dnew[index], intensity3Dnew[index], camNums.size(), ignoreCam);		

			// removing ambiguous and ghost particles
			Rem(pos3Dnew, intensity3Dnew, mindist_3D);

			// updating the reprojected image
			ReprojImage(pos3Dnew, OTFcalib, reproj, IPRflag);
//			ReprojImage(pos3Dnew, OTFcalib, reproj, 1.5);

			// updating the original image by removing correctly identified 3D particles
			for (int n = 0; n < camNums.size(); n++)
				for (int i = 0; i < Npixh; i++) {
					for (int j = 0; j < Npixw; j++) {
						int residual = (orig[camNums[n]][i][j] - config.fpt * reproj[camNums[n]][i][j]);
						orig[camNums[n]][i][j] = (residual < 0) ? 0 : residual;
					}
				}
		}
	} 
	filename.clear();
//	SaveParticlePositions(pos3Dnew.Get_PosDeque(), m_particle_position_addr + "frame" + to_string(frame) + "Loop" + to_string(loop_time) + ".txt");
	return pos3Dnew;
}



void IPR::ReprojImage(Frame matched3D, OTF& OTFcalib, deque<int**>& pixels_reproj, bool STB) {
	int size = psize;
	// doubling the area of reprojection for STB
	// also double for IPR for consistency with Shaking.
//	if (STB)
		size = 2 * psize;

	// intializing pixel_reproj to 0
	for (int camID = 0; camID < ncams; camID++) {
		for (int i = 0; i < Npixh; i++) {
			for (int j = 0; j < Npixw; j++) {
				pixels_reproj[camID][i][j] = 0;
			}
		}
	}

	Frame::const_iterator pIDend = matched3D.end();
	int p = 0;
	for (Frame::const_iterator pID = matched3D.begin(); pID != pIDend; ++pID) {

		// stores the particle's center on each camera image
		deque<Position> particle2Dcenters;

		for (int n = 0; n < ncams; n++) {
			// Use OTFparam matrix and interpolate (trilinear) to get the parameters at the current 3D particle location
			vector <double> otfParam = OTFcalib.OTFgrid(n, matched3D[p]); // otfParam contains a,b,c and alpha for camera 'n' at 3D position 'pos'
			// finding the 2D center
			Position pos2Dmm = camsAll[n].WorldToImage(*pID); pos2Dmm.Set_Z(0);
			particle2Dcenters.push_back(camsAll[n].Distort(pos2Dmm));

			// *Reporjecting* //
			// pixel range for each particle
			int xmin = max(1, (int)floor(particle2Dcenters[n].X() - size / 2));
			int ymin = max(1, (int)floor(particle2Dcenters[n].Y() - size / 2));
			int xmax = min(Npixw, (int)floor(particle2Dcenters[n].X() + size / 2));
			int ymax = min(Npixh, (int)floor(particle2Dcenters[n].Y() + size / 2));

			for (int x = xmin; x < xmax; x++) {
				for (int y = ymin; y < ymax; y++) {
					// reprojecting the particle using Gaussian ellipse
					int proj = round(PixelReproj(particle2Dcenters[n], otfParam, x, y));
					pixels_reproj[n][y][x] = max(pixels_reproj[n][y][x], proj);
					// important comment: Not sure max is the right thing to use here for overlapping particles
				}
			}
		}
		p++;
	}

}

void IPR::ReprojImage(Frame matched3D, OTF& OTFcalib, deque<int**>& pixels_reproj, double projsize) {
//	int size = psize;
//	// doubling the area of reprojection for STB
//	// also double for IPR for consistency with Shaking.
////	if (STB)
//		size = 2 * psize;
	projsize = projsize * psize;

	// intializing pixel_reproj to 0
	for (int camID = 0; camID < ncams; camID++) {
		for (int i = 0; i < Npixh; i++) {
			for (int j = 0; j < Npixw; j++) {
				pixels_reproj[camID][i][j] = 0;
			}
		}
	}

	Frame::const_iterator pIDend = matched3D.end();
	int p = 0;
	for (Frame::const_iterator pID = matched3D.begin(); pID != pIDend; ++pID) {

		// stores the particle's center on each camera image
		deque<Position> particle2Dcenters;

		for (int n = 0; n < ncams; n++) {
			// Use OTFparam matrix and interpolate (trilinear) to get the parameters at the current 3D particle location
			vector <double> otfParam = OTFcalib.OTFgrid(n, matched3D[p]); // otfParam contains a,b,c and alpha for camera 'n' at 3D position 'pos'
			// finding the 2D center
			Position pos2Dmm = camsAll[n].WorldToImage(*pID); pos2Dmm.Set_Z(0);
			particle2Dcenters.push_back(camsAll[n].Distort(pos2Dmm));

			// *Reporjecting* //
			// pixel range for each particle
			int xmin = max(1, (int)floor(particle2Dcenters[n].X() - projsize / 2));
			int ymin = max(1, (int)floor(particle2Dcenters[n].Y() - projsize / 2));
			int xmax = min(Npixw, (int)floor(particle2Dcenters[n].X() + projsize / 2));
			int ymax = min(Npixh, (int)floor(particle2Dcenters[n].Y() + projsize / 2));

			for (int x = xmin ; x < xmax; x++) {
				for (int y = ymin ; y < ymax; y++) {
					// reprojecting the particle using Gaussian ellipse
					int proj = round(PixelReproj(particle2Dcenters[n], otfParam, x, y));
					pixels_reproj[n][y][x] = max(pixels_reproj[n][y][x], proj);
					// important comment: Not sure max is the right thing to use here for overlapping particles
				}
			}
		}
		p++;
	}
}

// a function for Gaussian ellipse reprojection at position (x,y)
double IPR::PixelReproj(Position& particle2Dcenter, vector <double>& otfParam, int x, int y) {
	double xx = ((double)x - particle2Dcenter.X())*cos(otfParam[3]) + ((double)y - particle2Dcenter.Y())*sin(otfParam[3]);
	double yy = -((double)x - particle2Dcenter.X())*sin(otfParam[3]) + ((double)y - particle2Dcenter.Y())*cos(otfParam[3]);
	double value = otfParam[0] * exp(-(otfParam[1] * pow(xx, 2) + otfParam[2] * pow(yy, 2)));
	return(value);
}

// removing ghost and ambiguous particles
pair<int,int> IPR::Rem(Frame& pos3D, deque<double>& int3D, double mindist_3D) {
	int ambiguous = 0;

	// deleting particles that are very close to each other
	for (int i = 0; i < pos3D.NumParticles(); i++) {
		for (int j = i + 1; j < pos3D.NumParticles();) {
			if (Distance(pos3D[i], pos3D[j]) < 9 * mindist_3D * mindist_3D) {
				pos3D.Delete(j); int3D.erase(int3D.begin() + j);
				ambiguous++;
			}
			else
				j++;
		}
	}

	int ghost = 0;
	// deleting based on intensity
	double avgInt = 0;
	for (int i = 0; i < int3D.size(); i++) {
		if (int3D[i] > 0) {
			avgInt = avgInt + int3D[i];
		}
	}
	avgInt = avgInt / int3D.size();
//	cout<<"average intensity:"<<avgInt<<endl;

	for (int index = int3D.size() - 1; index >= 0; index--) {
		// removing a particle if its intensity falls below a % of the avg intensity
		if (int3D[index] < intensityLower * avgInt) {
			ghost3D.push_back(pos3D[index]);
			pos3D.Delete(index); int3D.erase(int3D.begin() + index);
			ghost++;
		}
	}

	return make_pair(ambiguous, ghost);
}

void IPR::FullData(Position& pos, double intensity, int Cams, int ignoreCam) {
	deque< deque<double> > pos2D(4);
	for (int camID = 0; camID < 4; camID++) {
		if (camID < Cams) {
			deque<double> tmp2D(2);
			Position tmp = camsAll[camID].WorldToImage(pos); tmp.Set_Z(0);
			tmp = camsAll[camID].Distort(tmp);
			tmp2D[0] = tmp.X(); tmp2D[1] = tmp.Y();
			pos2D[camID] = tmp2D;
		}
		else {
			/*
			 * Modified by Shiyong Tan, 1/31/18
			 * Illegal initialization with {...}
			 * Start:
			 */
//			deque<double> tmp = { 0.0, 0.0 };
			deque<double> tmp(2);
			tmp[0] = 0.0; tmp[1] = 0.0;
			// End
			pos2D[camID] = tmp;
		}
	}
	Position temp(pos.X(), pos.Y(), pos.Z(), pos2D[0][0], pos2D[0][1], pos2D[1][0], pos2D[1][1], pos2D[2][0], pos2D[2][1], pos2D[3][0], pos2D[3][1], intensity);
	pos = temp;
}

// creates a matfile of positions
void IPR::SaveParticlePositions(deque<Position> pos, string file_path) {
/*
 * Modified by Shiyong Tan, 2/6/18
 * Discard using matio, use DataIO instead
 * Start:
 */
//	// Create a .mat file with pos3D
//	size_t sizeofpos3D = pos.size();
//	mat_t    *matfp;
//	matvar_t *cell_array, *cell_element;
//	size_t dims[2] = { sizeofpos3D, 1 };
//	stringstream s; s << "pos3Dframe" << frame;
//	string mat_name = name + ".mat";
//	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
//	switch (NULL == matfp) {
//		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
//		break;
//	}
//
//	cell_array = Mat_VarCreate(s.str().c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//	if (NULL == cell_array) {
//		fprintf(stderr, "Error creating variable for 'pos3D'\n");
//	}
//	else {
//		for (int i = 0; i < sizeofpos3D; i++) {
//			dims[0] = 1;
//			dims[1] = 12;
//			double temp[12] = { pos[i].X(),pos[i].Y(),pos[i].Z(),pos[i].X1(),pos[i].Y1(),pos[i].X2(),pos[i].Y2(),pos[i].X3(),pos[i].Y3(),pos[i].X4(),pos[i].Y4(),pos[i].Info() };
//			cell_element = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, temp, 0);
//			switch (NULL == cell_element) {
//				fprintf(stderr, "Error creating cell element variable\n");
//				Mat_VarFree(cell_array);
//				Mat_Close(matfp);
//				break;
//			}
//			Mat_VarSetCell(cell_array, i, cell_element);
//		}
//	}
//
//	Mat_VarWrite(matfp, cell_array, MAT_COMPRESSION_NONE);
//	Mat_VarFree(cell_array);
//	Mat_Close(matfp);

	// TODO: Check whether it works. Shiyong Tan
	size_t sizeofpos3D = pos.size();
	double* array = new double[sizeofpos3D * 12];
	//Convert Position into array
	Position2Array(pos, array);

	NumDataIO<double> data_io;
	data_io.SetFilePath(file_path);  // setting the path to save the data.
	data_io.SetTotalNumber(sizeofpos3D * 12);
	data_io.WriteData((double*) array);
	delete[] array;
// END
}

Frame IPR::ReadParticlePositions(string file_path) {
	NumDataIO<double> data_io;
	data_io.SetFilePath(file_path);
	int num = data_io.GetTotalNumber(); //Get total data number
	double array[num / 12][12];
	data_io.ReadData((double*) array);
	deque<Position> pos = Array2Position(num / 12, array);
	Frame frame(pos);
	return frame;
}

// creates a matfile for images
void IPR::MatfileImage(deque<int**>& pix, string name) {

	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Matio is discared. Use DataIO instead.
	 * Start:
	 */
//	size_t sizeofpos3D = pix.size();
//	mat_t    *matfp;
//	size_t dims[2] = { Npixw, 1 };
//	string mat_name = name + ".mat";
//	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
//	switch (NULL == matfp) {
//		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
//		break;
//	}
//
//	for (int i = 0; i < sizeofpos3D; i++) {
//		dims[0] = Npixw;
//		dims[1] = 1;
//		stringstream cam;
//		cam << "cam" << i;
//		matvar_t **cell_element = new matvar_t*[Npixw];
//		matvar_t *cell_array = Mat_VarCreate(cam.str().c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//		if (NULL == cell_array) {
//			fprintf(stderr, "Error creating variable for 'pos3D'\n");
//		}
//		else {
//			for (int j = 0; j < Npixw; j++) {
//				dims[0] = 1;
//				dims[1] = Npixh;
//				int* temp = new int[Npixh];
//				temp = pix[i][j];
//
//				cell_element[j] = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims, temp, 0);
//
//				switch (NULL == cell_element) {
//					fprintf(stderr, "Error creating cell element variable\n");
//					Mat_VarFree(cell_array);
//					Mat_Close(matfp);
//					break;
//				}
//				Mat_VarSetCell(cell_array, j, cell_element[j]);
//				//Mat_VarFree(cell_element);
//				//delete[] temp;
//			}
//		}
//		Mat_VarWrite(matfp, cell_array, MAT_COMPRESSION_NONE);
//		for (int j = 0; j < Npixw; j++) {
//			Mat_VarFree(cell_element[j]);
//		}
//		delete[] cell_element;
//		//Mat_VarFree(cell_array);
//	}
//	Mat_Close(matfp);
	// TODO: check whether it works.
	size_t sizeofpos3D = pix.size();
	// get all the values of pixels into 3D matrix
	int tmp[sizeofpos3D][Npixw][Npixh];
	for (int i = 0; i < sizeofpos3D; i++) {
		for (int j = 0; j < Npixw; j++ ) {
			for (int k = 0; k < Npixh; k++) tmp[i][j][k] = pix[i][j][k];
		}
	}
	NumDataIO<int> data_io;
	data_io.SetFilePath(name + ".txt");
	data_io.SetTotalNumber(sizeofpos3D * Npixw * Npixh);
	data_io.WriteData((int*) tmp);
// End
}

void IPR::Load_2Dpoints(string path, int frame, int ignoreCam) {
	/*
	 * Modified by Shiyong Tan, 2/7/18
	 * Discard using matio, use DataIo istead.
	 * Start:
	 */
//	stringstream s; s << path << ".mat";
//	string file = s.str();
//
//	const char *fileName = file.c_str();
//	mat_t *mat = Mat_Open(fileName, MAT_ACC_RDONLY);
//
//	if (mat == NULL) {
//		cout << " error in reading the 2D pos mat file" << endl;
//	}
//
//	for (int id = 0; id < ncams; id++) {
//		if (id != ignoreCam) {
//			int k = id;
//
//			//if (id == 2)
//			//	k = 3;
//			//if (id == 3)
//			//	k = 2;
//			stringstream s; s << "frame_" << frame;
//
//			string varName = s.str();
//			Frame pos2D;
//
//			if (mat) {
//				//std::cout << "Open file to read\n\tmat == " << mat << "\n";
//
//				matvar_t *matVar = 0;
//				matVar = Mat_VarRead(mat, (char*)varName.c_str());
//
//				if (matVar) {
//					int rows;	int cols;
//					rows = matVar->dims[0]; cols = matVar->dims[1];
//					unsigned namesize = matVar->nbytes / matVar->data_size;
//					double *namedata = static_cast<double*>(matVar->data);
//					for (int i = 0; i < rows; i++) {
//						Position pos(namedata[(2*id)*rows + i], namedata[(2*id + 1)*rows + i], 0);
//						pos2D.Add(pos);
//					}
//				}
//				else
//					cout << "Cannot open mat variable\n";
//			}
//			else
//				cout << "Cannot open mat file\n";
//
//			iframes.push_back(pos2D);
//		}
//	}
//
//	Mat_Close(mat);
	//TODO: to check whether it works.
	NumDataIO<double> data_io;
	data_io.SetFilePath(path + ".txt");
	int total_number = data_io.GetTotalNumber(); //get the total number of elements in txt file
	// suppose the data format is: X1, Y1, X2, Y2, X3, Y3, X4, Y4
	int rows = total_number / 8;
	double points_array[rows][8];
	data_io.ReadData((double*) points_array);
	for (int id = 0; id < ncams; id++) {
		Frame pos2D;
		if (id != ignoreCam) {
			for (int i = 0; i < rows; i++) {
				Position pos(points_array[i][2 * id], points_array[i][2 * id + 1], 0);
				pos2D.Add(pos);
			}
		}
		iframes.push_back(pos2D);
	}
	// End
}

void IPR::Position2Array(deque<Position> pos, double* array) {
	size_t sizeofpos3D = pos.size();
	for (int i = 0; i < sizeofpos3D; i++) {
		array[i * 12 + 0] = pos[i].X(); array[i * 12 + 1] =pos[i].Y(); array[i * 12 + 2] = pos[i].Z();
		array[i * 12 + 3] = pos[i].X1(); array[i * 12 + 4] = pos[i].Y1(); array[i * 12 + 5] = pos[i].X2(); array[i * 12 + 6] =pos[i].Y2();
		array[i * 12 + 7] = pos[i].X3(); array[i * 12 + 8] = pos[i].Y3(); array[i * 12 + 9] = pos[i].X4(); array[i * 12 + 10] = pos[i].Y4();
		array[i * 12 + 11] = pos[i].Info();
	}

}

deque<Position> IPR::Array2Position(int num_particle, double array[][12]) {
	deque<Position> pos;
	for(int i = 0; i < num_particle; i++) {
		Position point(array[i][0], array[i][1], array[i][2], array[i][3], array[i][4], array[i][5],  array[i][6],
				array[i][7], array[i][8], array[i][9], array[i][10], array[i][11]);
		pos.push_back(point);
	}
	return pos;
}

//Save Particle intensity
void IPR::SaveParticleIntensity(deque<double> intensity, string file_path) {
	NumDataIO<double> data_io;
	int num = intensity.size();
	data_io.SetFilePath(file_path);
	double intensity_array[num]; // Convert deque into double
	for (int i = 0; i <  num; i++) {
		intensity_array[i] = intensity[i];
	}
	data_io.SetTotalNumber(num);
	data_io.WriteData((double*) intensity_array);
}

// Read the particle intensity
deque<double> IPR::ReadParticleIntensity(string file_path) {
	NumDataIO<double> data_io;
	data_io.SetFilePath(file_path);
	int num = data_io.GetTotalNumber();
	double intensity_array[num];
	data_io.ReadData((double*) intensity_array);
	deque<double> intensity;
	for (int i = 0; i < num; i++) {
		intensity.push_back(intensity_array[i]);
	}
	return intensity;
}

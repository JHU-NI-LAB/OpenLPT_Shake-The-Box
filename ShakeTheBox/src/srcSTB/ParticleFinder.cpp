/*
 *  ParticleFinder.cpp
 *
 *  This is the implementation file for ParticleFinder objects.
 *
 *  Last update: 10/12/11 by NTO
 *
 */

#include <iostream>
#include <fstream>
#include <queue>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include<sstream>
#include<string>

#include <ParticleFinder.h>
#include <Logs.h>
#include <Position.h>
#include "Common.h"

/*
 * Modified by Shiyong Tan, 2/5/18
 * The matio library is discarded and use DataIO instead.
 * Start:
 */
//#include <matio.h>
#include "NumDataIO.h"
// End

using namespace std;

//ParticleFinder::ParticleFinder(int**& p, int rows, int cols, int depth, int threshold)
//: pixels(p), rows(rows), cols(cols)
//{
//	int colors = depth;
//
//  // walk through the given array, skipping the first and last row and column
//  for (int i = 1; i < (rows - 1); ++i) {
//    for (int j = 1; j < (cols - 1); ++j) {
//      // is this pixel a local maximum above threshold?
//      if ((pixels[i][j] >= threshold) && IsLocalMax(i, j)) {
//				// read in the local maximum pixel value as well as the values to its
//				// left and right and top and bottom in order to calculate the
//				// particle center.  note: add 0.5 to row and column values to put
//				// the pixel origin in its center
//		  double x1 = (j - 1);// +0.5;
//		  double x2 = j;// +0.5;
//		  double x3 = (j + 1);// +0.5;
//		  double y1 = (i - 1);// +0.5;
//		  double y2 = i;// + 0.5;
//		  double y3 = (i + 1);// +0.5;
//
//				// check the pixel values to make sure we have no corrupted memory
//				if ((abs(pixels[i][j]) > colors) ||
//						(abs(pixels[i - 1][j]) > colors) ||
//						(abs(pixels[i + 1][j]) > colors) ||
//						(abs(pixels[i][j - 1]) > colors) ||
//						(abs(pixels[i][j + 1]) > colors)) {
//					// this is a serious problem: we can't continue
//					stringstream s;
//					s << "Pixel [" << i << "," << j << "] out of range with " << colors;
//					/*
//					 * Modified by Shiyong Tan, 2/5/18
//					 * Update the IO interface.
//					 * Start:
//					 */
//					MatfileImage(pixels, "ErrorImage");
//					// End
//					throw out_of_range(s.str());
//				}
//
//				// find the column value, moving 0 intensities to 0.0001 so we can
//				// take their log
//				double lnz1, lnz2, lnz3;
//				if (colors == 255) {
//					if (pixels[i][j - 1] == 0) {
//						lnz1 = log(0.0001);
//					} else {
//						lnz1 = Logs::log8bit[pixels[i][j - 1]];
//					}
//					if (pixels[i][j] == 0) {
//            			lnz2 = log(0.0001);
//					} else {
//						lnz2 = Logs::log8bit[pixels[i][j]];
//					}
//					if (pixels[i][j + 1] == 0) {
//						lnz3 = log(0.0001);
//					} else {
//						lnz3 = Logs::log8bit[pixels[i][j + 1]];
//					}
//        } else if (colors == 65535) {
//          lnz1 = Logs::log16bit[pixels[i][j - 1]];
//          lnz2 = Logs::log16bit[pixels[i][j]];
//					lnz3 = Logs::log16bit[pixels[i][j + 1]];
//				} else {
//          if (pixels[i][j - 1] == 0) {
//						lnz1 = log(0.0001);
//					} else {
//						lnz1 = log(static_cast<double>(pixels[i][j - 1]));
//					}
//					if (pixels[i][j] == 0) {
//            lnz2 = log(0.0001);
//					} else {
//						lnz2 = log(static_cast<double>(pixels[i][j]));
//					}
//					if (pixels[i][j + 1] == 0) {
//						lnz3 = log(0.0001);
//					} else {
//						lnz3 = log(static_cast<double>(pixels[i][j + 1]));
//					}
//        }
//
//				double xc = -0.5 * ((lnz1 * ((x2 * x2) - (x3 * x3))) - (lnz2 * ((x1 * x1) - (x3 * x3))) + (lnz3 * ((x1 * x1) - (x2 * x2)))) / ((lnz1 * (x3 - x2)) - (lnz3 * (x1 - x2)) + (lnz2 * (x1 - x3)));
//
//				// were these numbers valid?
//				if (!isfinite(xc)) {
//          // no -- we had a problem.  drop this particle
//					continue;
//				}
//
//				// find the row value
//				if (colors == 255) {
//          if (pixels[i - 1][j] == 0) {
//						lnz1 = log(0.0001);
//					} else {
//						lnz1 = Logs::log8bit[pixels[i - 1][j]];
//					}
//					if (pixels[i + 1][j] == 0) {
//            lnz3 = log(0.0001);
//					} else {
//						lnz3 = Logs::log8bit[pixels[i + 1][j]];
//					}
//				} else if (colors == 65535) {
//					lnz1 = Logs::log16bit[pixels[i - 1][j]];
//					lnz3 = Logs::log16bit[pixels[i + 1][j]];
//				} else {
//          if (pixels[i - 1][j] == 0) {
//						lnz1 = log(0.0001);
//					} else {
//						lnz1 = log(static_cast<double>(pixels[i - 1][j]));
//					}
//					if (pixels[i + 1][j] == 0) {
//            lnz3 = log(0.0001);
//					} else {
//						lnz3 = log(static_cast<double>(pixels[i + 1][j]));
//					}
//				}
//
//				double yc = -0.5 * ((lnz1 * ((y2 * y2) - (y3 * y3))) - (lnz2 * ((y1 * y1) - (y3 * y3))) +
//														(lnz3 * ((y1 * y1) - (y2 * y2))))
//														/ ((lnz1 * (y3 - y2)) - (lnz3 * (y1 - y2)) + (lnz2 * (y1 - y3)));
//
//				// check these numbers too
//				if (!isfinite(yc)) {
//					// a problem occurred, so we'll drop these numbers and move on
//					continue;
//				}
//				//cout << "xc: " << xc << " yc: " << yc << endl;
//			y.push_back(yc);
//			x.push_back(xc);
//      }
//    }
//  }
//}

void ParticleFinder::GetParticle2DCenter(int depth, int threshold) {
	int colors = depth;

  // walk through the given array, skipping the first and last row and column
  for (int i = 1; i < (rows - 1); ++i) {
    for (int j = 1; j < (cols - 1); ++j) {
      // is this pixel a local maximum above threshold?
      if ((pixels[i][j] >= threshold) && IsLocalMax(i, j)) {
				// read in the local maximum pixel value as well as the values to its
				// left and right and top and bottom in order to calculate the
				// particle center.  note: add 0.5 to row and column values to put
				// the pixel origin in its center
		  double x1 = (j - 1);// +0.5;
		  double x2 = j;// +0.5;
		  double x3 = (j + 1);// +0.5;
		  double y1 = (i - 1);// +0.5;
		  double y2 = i;// + 0.5;
		  double y3 = (i + 1);// +0.5;

				// check the pixel values to make sure we have no corrupted memory
				if ((abs(pixels[i][j]) > colors) ||
						(abs(pixels[i - 1][j]) > colors) ||
						(abs(pixels[i + 1][j]) > colors) ||
						(abs(pixels[i][j - 1]) > colors) ||
						(abs(pixels[i][j + 1]) > colors)) {
					// this is a serious problem: we can't continue
					stringstream s;
					s << "Pixel [" << i << "," << j << "] out of range with " << colors;
					/*
					 * Modified by Shiyong Tan, 2/5/18
					 * Update the IO interface.
					 * Start:
					 */
					SaveImage(pixels, "ErrorImage");
					// End
					throw out_of_range(s.str());
				}

				// find the column value, moving 0 intensities to 0.0001 so we can
				// take their log
				double lnz1, lnz2, lnz3;
				if (colors == 255) {
					if (pixels[i][j - 1] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = Logs::log8bit[pixels[i][j - 1]];
					}
					if (pixels[i][j] == 0) {
            			lnz2 = log(0.0001);
					} else {
						lnz2 = Logs::log8bit[pixels[i][j]];
					}
					if (pixels[i][j + 1] == 0) {
						lnz3 = log(0.0001);
					} else {
						lnz3 = Logs::log8bit[pixels[i][j + 1]];
					}
        } else if (colors == 65535) {
          lnz1 = Logs::log16bit[pixels[i][j - 1]];
          lnz2 = Logs::log16bit[pixels[i][j]];
					lnz3 = Logs::log16bit[pixels[i][j + 1]];
				} else {
          if (pixels[i][j - 1] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = log(static_cast<double>(pixels[i][j - 1]));
					}
					if (pixels[i][j] == 0) {
            lnz2 = log(0.0001);
					} else {
						lnz2 = log(static_cast<double>(pixels[i][j]));
					}
					if (pixels[i][j + 1] == 0) {
						lnz3 = log(0.0001);
					} else {
						lnz3 = log(static_cast<double>(pixels[i][j + 1]));
					}
        }

				double xc = -0.5 * ((lnz1 * ((x2 * x2) - (x3 * x3))) - (lnz2 * ((x1 * x1) - (x3 * x3))) + (lnz3 * ((x1 * x1) - (x2 * x2)))) / ((lnz1 * (x3 - x2)) - (lnz3 * (x1 - x2)) + (lnz2 * (x1 - x3)));

				// were these numbers valid?
				if (!isfinite(xc)) {
          // no -- we had a problem.  drop this particle
					continue;
				}

				// find the row value
				if (colors == 255) {
          if (pixels[i - 1][j] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = Logs::log8bit[pixels[i - 1][j]];
					}
					if (pixels[i + 1][j] == 0) {
            lnz3 = log(0.0001);
					} else {
						lnz3 = Logs::log8bit[pixels[i + 1][j]];
					}
				} else if (colors == 65535) {
					lnz1 = Logs::log16bit[pixels[i - 1][j]];
					lnz3 = Logs::log16bit[pixels[i + 1][j]];
				} else {
          if (pixels[i - 1][j] == 0) {
						lnz1 = log(0.0001);
					} else {
						lnz1 = log(static_cast<double>(pixels[i - 1][j]));
					}
					if (pixels[i + 1][j] == 0) {
            lnz3 = log(0.0001);
					} else {
						lnz3 = log(static_cast<double>(pixels[i + 1][j]));
					}
				}

				double yc = -0.5 * ((lnz1 * ((y2 * y2) - (y3 * y3))) - (lnz2 * ((y1 * y1) - (y3 * y3))) +
														(lnz3 * ((y1 * y1) - (y2 * y2))))
														/ ((lnz1 * (y3 - y2)) - (lnz3 * (y1 - y2)) + (lnz2 * (y1 - y3)));

				// check these numbers too
				if (!isfinite(yc)) {
					// a problem occurred, so we'll drop these numbers and move on
					continue;
				}
				//cout << "xc: " << xc << " yc: " << yc << endl;
			y.push_back(yc);
			x.push_back(xc);
      }
    }
  }
}

void ParticleFinder::WriteToFile(string filename) {
  // now we have all the particle centers and can write them to a file
  ofstream outfile(filename.c_str(), ios::out);
  outfile << "# Modified Gaussian fitting (3-point method)\n";
	
  deque<double>::const_iterator xi_end = x.end();
  deque<double>::const_iterator yi = y.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
    outfile << *xi << "\t" << *yi << "\n";
  }
	
  outfile.close();
}

Frame ParticleFinder::CreateFrame() {
	deque<Position> pos;
	deque<double>::const_iterator xi_end = x.end();
	deque<double>::const_iterator yi = y.begin();
	for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
		pos.push_back(Position(*xi, *yi, 0));
	}
	return Frame(pos);
}

void ParticleFinder::Squash(double rad) {
  if (rad < 1) {
    // don't do anything for small cluster radii, for efficiency
    return;
  }
  
  deque<deque<double>::const_iterator> bad;
  deque<double> newx;
  deque<double> newy;
  
  deque<double>::const_iterator xi_end = x.end();
  deque<double>::const_iterator yi = y.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
    // have we looked at this entry before?
    if (find(bad.begin(), bad.end(), xi) != bad.end()) {
      // yes: skip it.
      continue;
    }
    
    double avgx = *xi;
    double avgy = *yi;
    int N = 1;
    // loop over other positions, looking for neighbors
    deque<double>::const_iterator xi2 = xi + 1;
    deque<double>::const_iterator yi2 = yi + 1;
    for (; xi2 != xi_end; ++xi2, ++yi2) {
      if (pow(*xi-*xi2,2)+pow(*yi-*yi2,2) > rad*rad) {
        continue;
      }
      // these positions are nearby
      avgx += *xi2;
      avgy += *yi2;
      ++N;
    }
    // did we find a cluster?
    if (N == 1) {
      // nope.
      continue;
    }
    
    bad.push_back(xi);
    int oldN = N;
    
    // now loop through again, looking for everything near *this average*
    // location, and re-compute the average. continue until we converge.
    while (true) {
      double ax = avgx / static_cast<double>(N);
      double ay = avgy / static_cast<double>(N);
      avgx = *xi;
      avgy = *yi;
      N = 1;
      deque<double>::const_iterator xi2 = xi + 1;
      deque<double>::const_iterator yi2 = yi + 1;
      for (; xi2 != xi_end; ++xi2, ++yi2) {
        if (pow(ax-*xi2,2)+pow(ay-*yi2,2) > rad*rad) {
          continue;
        }
        // these positions are nearby
        avgx += *xi2;
        avgy += *yi2;
        ++N;
      }
      if (oldN == N) {
        // we converged! now do bookkeeping
        newx.push_back(avgx / static_cast<double>(N));
        newy.push_back(avgy / static_cast<double>(N));
        deque<double>::const_iterator xi2 = xi + 1;
        deque<double>::const_iterator yi2 = yi + 1;
        for (; xi2 != xi_end; ++xi2, ++yi2) {
          if (pow(ax-*xi2,2)+pow(ay-*yi2,2) <= rad*rad) {
            bad.push_back(xi2);
          }
        }
        break;
      } else {
        // try again
        oldN = N;
      }
    }
  }
  
  // finally, find the good positions left in the original queue
  yi = y.begin();
  for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
    if (find(bad.begin(), bad.end(), xi) == bad.end()) {
      newx.push_back(*xi);
      newy.push_back(*yi);
      //cout << *xi << "\t" << *yi << endl;
    }
  }
  x = newx;
  y = newy;
}

bool ParticleFinder::IsLocalMax(int r, int c)
{
  int val = pixels[r][c];

  // check only the 4 neighbors directly above, below, to the left, and to the right
  int num_pix = 4;
//  if (pixels[r][c - 1] > val) --num_pix;
//  if (pixels[r][c + 1] > val) --num_pix;
//  if (pixels[r - 1][c] > val) --num_pix;
//  if (pixels[r + 1][c] > val) --num_pix;
//
//  if (num_pix <= 2) return false;

  if ((pixels[r][c - 1] > val) || (pixels[r][c + 1] > val) || (pixels[r - 1][c] > val) || (pixels[r + 1][c] > val)) {
    return false;
  }
	
  return true;
}


/*
 * Modified by Shiyong Tan, 2/5/18
 * The matio library is discarded and use DataIO instead.
 * Start:
 */
// creates a matfile for images
//void ParticleFinder::MatfileImage(int** pix, string name) {
//
//	size_t sizeofpos3D = 1;
//	mat_t    *matfp;
//	size_t dims[2] = { rows, 1 };
//	string mat_name = name + ".mat";
//	matfp = Mat_CreateVer(mat_name.c_str(), NULL, MAT_FT_DEFAULT);
//	switch (NULL == matfp) {
//		fprintf(stderr, "Error creating MAT file \"pos3D.mat\"!\n");
//		break;
//	}
//
//	for (int i = 0; i < sizeofpos3D; i++) {
//		dims[0] = rows;
//		dims[1] = 1;
//		matvar_t *cell_array, *cell_element;
//		cell_array = Mat_VarCreate("Image", MAT_C_CELL, MAT_T_CELL, 2, dims, NULL, 0);
//		if (NULL == cell_array) {
//			fprintf(stderr, "Error creating variable for 'pos3D'\n");
//		}
//		else {
//			for (int j = 0; j < rows; j++) {
//				dims[0] = 1;
//				dims[1] = rows;
//				int* temp = new int[cols];
//				temp = pix[j];
//
//				cell_element = Mat_VarCreate(NULL, MAT_C_INT32, MAT_T_INT32, 2, dims, temp, 0);
//
//				switch (NULL == cell_element) {
//					fprintf(stderr, "Error creating cell element variable\n");
//					Mat_VarFree(cell_array);
//					Mat_Close(matfp);
//					break;
//				}
//				Mat_VarSetCell(cell_array, j, cell_element);
//
//			}
//		}
//		Mat_VarWrite(matfp, cell_array, MAT_COMPRESSION_NONE);
//		//		Mat_VarFree(cell_array);
//	}
//	Mat_Close(matfp);
//}

//TODO check whether it works Shiyong Tan 2/5/18
void ParticleFinder::SaveImage(int** pix, string name) {
	NumDataIO<int> data_io;
	data_io.SetFilePath(name + ".txt");
	data_io.SetTotalNumber(rows * cols);
	data_io.WriteData((int*) pix);
}
//End

// Save the particle 2D center to txt file
void ParticleFinder::SaveParticle2DCenter(string file_path) {
	NumDataIO<double> data_io;
	file_path.erase(file_path.end() - 3, file_path.end()); // erase "tif"
	data_io.SetFilePath(file_path + "txt");
	deque<double>::const_iterator xi_end = x.end();
	deque<double>::const_iterator yi = y.begin();
	double particle_2Dcenter[x.size()][2];
	int i = 0;
	for (deque<double>::const_iterator xi = x.begin(); xi != xi_end; ++xi, ++yi) {
		particle_2Dcenter[i][0] = *xi;
		particle_2Dcenter[i][1] = *yi;
		i++;
	}
	data_io.SetTotalNumber(i * 2);
	data_io.WriteData((double *) particle_2Dcenter);
}

// Read the particle 2D center to txt file
Frame ParticleFinder::ReadParticle2DCenter(string file_path) {
	NumDataIO<double> data_io;
	file_path.erase(file_path.end() - 3, file_path.end()); // erase "tif"
	data_io.SetFilePath(file_path + "txt");
	int num = data_io.GetTotalNumber();
	double* particle_2Dcenter = new double[num];
	deque<Position> pos;
	data_io.SetTotalNumber(num);
	data_io.ReadData(particle_2Dcenter);
	for (int i = 0; i < num / 2; i++) {
		pos.push_back(Position(particle_2Dcenter[i*2], particle_2Dcenter[i*2 + 1], 0));
	}
	return pos;
}


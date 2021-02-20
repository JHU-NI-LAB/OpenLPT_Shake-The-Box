/*
 *  ParticleFinder.h
 * 
 *  A ParticleFinder object knows how to find particles given a pixel array.
 *  Currently, it does this by using the 1D Gaussian Estimator procedure 
 *  described in:
 *
 *  N.T. Ouellette, H. Xu, and E. Bodenschatz, "A quantitative study of three-
 *  dimensional Lagrangian particle tracking algorithms," Exp. Fluids 40, 301-313 (2006).
 *
 *  To change the particle-finding method, re-implement the constructor.
 *
 *  Last update: 10/12/11 by NTO
 *
 */

#ifndef PARTICLEFINDER_H
#define PARTICLEFINDER_H

#include <string>
#include <deque>
#include <stdexcept>

#include <Frame.h>

class ParticleFinder {

public:
  // constructor: process the given pixel array
  ParticleFinder(int**& p, int rows, int cols) : pixels(p), rows(rows), cols(cols) {};
  // destructor
  ~ParticleFinder() {};

  /*
   * Function: to get the 2D center of particles in a image
   * Input: depth:
   * 		threshold:
   * Output: None.
   */
  void GetParticle2DCenter(int depth, int threshold);

  // write the particle centers to a file
  void WriteToFile(std::string filename);
  // make a Frame object with the particle positions.
  Frame CreateFrame();
   
  // return the number of particles found
  int NumParticles() const;
  
  // remove large clusters
  void Squash(double rad);

  void SaveImage(int** pix, std::string name);

  /*
   * Function: to save the 2D center of particles in a image to txt file
   * Input: file_path: the path to save the data
   * Output: None.
   */
  void SaveParticle2DCenter(std::string file_path);

  /*
   * Function: to read the 2D center of particles in a image from txt file
   * Input: file_path: the path to read the data
   * Output: a frame format data with a list of particle 2D centers
   */
  Frame ReadParticle2DCenter(std::string file_path);

private:
	int** pixels;
	int rows;
	int cols;
	
  // store vectors of the x and y coordinates of the particles
  std::deque<double> x;
  std::deque<double> y;

	// helper function
  bool IsLocalMax(int r, int c);

};

inline int ParticleFinder::NumParticles() const
{
  return x.size();
}

#endif // PARTICLEFINDER_H

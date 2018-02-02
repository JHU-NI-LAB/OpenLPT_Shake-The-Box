#ifndef Tiff2DFinder_H
#define	Tiff2DFinder_H

#include <string>
#include <deque>
#include <stdexcept>
#include <utility>

#include <Frame.h>
#include <Image.h>
#include <Tiffload.h>
#include <ParticleFinder.h>

using namespace std;

class Tiff2DFinder {

public:
	// constructor: takes # of cams, threshold and tif filename for each cam
	Tiff2DFinder(int ncams, int thresh, deque<string>& f); //throw(runtime_error);
	// destructor
	~Tiff2DFinder() {};

	//template <class my_int>
	void Particle2DList(deque<int**>& pixels, deque<Frame>& iframes) throw(runtime_error);
	void FillPixels(deque<int**>& pixels) throw(runtime_error);

	// Get colors
	int Get_colors();
private:

	deque<string> filename;
	int ncams;
	int threshold;
	int Npixh;
	int Npixw;
	int colors;
};

inline int Tiff2DFinder::Get_colors() {
	return colors;
}



#endif

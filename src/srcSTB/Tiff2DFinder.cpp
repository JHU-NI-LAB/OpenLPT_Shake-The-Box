#include <cstdlib>
#include <string>

#include <Position.h>
#include <Frame.h>
#include <Image.h>
#include <Tiffload.h>
#include <Tiff2DFinder.h>
#include <ParticleFinder.h>

using namespace std;

Tiff2DFinder::Tiff2DFinder(int cams, int thresh, deque<string>& f) :
	ncams(cams), threshold(thresh), filename(f)
{
	if (filename.size() != ncams) {
		throw runtime_error("Number of Tiff names should be equal to number of cams in Tiff2DFinder");
	}
}

//template <class my_int>
// a function to fill pixels with the px from actual images and fill iframes with 2D centers from each particle
void Tiff2DFinder::Particle2DList(deque<int**>& pixels, deque<Frame>& iframes) throw(runtime_error)
{
	deque<Image*> images;
	deque<Image*>::iterator cur;
	for (int i = 0; i < filename.size(); ++i) {
		images.push_back(new Tiffload(filename[i]));
	}
	
	int camID = 0;
		
	for (cur = images.begin(); cur != images.end(); ++cur)
	{
		try {
			(*cur)->Open();
		}
		catch (invalid_argument& e) {
			cerr << e.what() << endl;
			throw runtime_error("Caught invalid_argument in Tiff2DFinder::Particle2DList()");
		}
		catch (runtime_error& e) {
			cerr << e.what() << endl;
			throw runtime_error("Caught runtime_error in Tiff2DFinder::Particle2DList()");
		}

		Npixh = (*cur)->GetRows();
		Npixw = (*cur)->GetCols();
		colors = (*cur)->GetDepth();

		try {
			pixels.push_back(new int*[Npixh]);
			for (int i = 0; i < Npixh; ++i) {
				pixels[camID][i] = new int[Npixw];
			}
		}
		catch (bad_alloc&) {
			throw runtime_error("Failure to allocate pixels in ImageSequence::Particle2DList()");
		}

		try {
			(*cur)->GetPixels(pixels[camID]);
		}
		catch (runtime_error& r) {
			cerr << r.what() << endl;
			throw runtime_error("Error in ImageSequence::Particle2DList()");
		}

		try {
			ParticleFinder p(pixels[camID], Npixh, Npixw);//, colors, threshold);
			p.GetParticle2DCenter(colors, threshold);
			iframes.push_back(p.CreateFrame());
		}
		catch (out_of_range& e) {
			cerr << e.what() << endl;
			throw runtime_error("Caught out_of_range in ImageSequence::Particle2DList()");
		}
		camID++;
	}

	for (int i = 0; i < filename.size(); ++i) {
		delete[] images[i];
	}

}

//  function to fill pixels with px from actual images
void Tiff2DFinder::FillPixels(deque<int**>& pixels) throw(runtime_error)
{
	deque<Image*> images;
	deque<Image*>::iterator cur;
	for (int i = 0; i < filename.size(); ++i) {
		images.push_back(new Tiffload(filename[i]));
	}

	int camID = 0;

	for (cur = images.begin(); cur != images.end(); ++cur)
	{
		try {
			(*cur)->Open();
		}
		catch (invalid_argument& e) {
			cerr << e.what() << endl;
			throw runtime_error("Caught invalid_argument in Tiff2DFinder::Particle2DList()");
		}
		catch (runtime_error& e) {
			cerr << e.what() << endl;
			throw runtime_error("Caught runtime_error in Tiff2DFinder::Particle2DList()");
		}

		Npixh = (*cur)->GetRows();
		Npixw = (*cur)->GetCols();
		colors = (*cur)->GetDepth();

		try {
			(*cur)->GetPixels(pixels[camID]);
		}
		catch (runtime_error& r) {
			cerr << r.what() << endl;
			throw runtime_error("Error in ImageSequence::Particle2DList()");
		}

		camID++;
	}

	for (int i = 0; i < filename.size(); ++i) {
		delete images[i];  //TODO: comments for debug.
//		images.erase(images.begin() + i);
	}
}

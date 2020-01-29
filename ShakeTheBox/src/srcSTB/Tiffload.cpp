/* 
 *  Tiff.cpp
 *
 *  Implementation file for Tiff objects.
 *
 *  Last update: 10/2/09 by NTO
 *
 */

// tiffio.h is provided by the libtiff package
#include <cstdlib>
#include <iostream>

#include <tiffio.h>
//#include <tiff.h>
#include <Tiffload.h>

using namespace std;

Tiffload::Tiffload(string f) throw(invalid_argument, runtime_error) 
	: filename(f)
{
	try {
	  Open();
	} catch (invalid_argument& e) {
		cerr << e.what() << endl;
		throw invalid_argument("Caught invalid_argument in Tiff constructor");
	} catch (runtime_error& e) {
		cerr << e.what() << endl;
		throw runtime_error("Caught runtime_error in Tiff constructor");
	}
}

void Tiffload::Open() throw(invalid_argument, runtime_error)
{
	// have we already opened the image?
	if (buffer != NULL) {
		return;
	}
	
  TIFF* image;
  if ((image = TIFFOpen(filename.c_str(), "r")) == NULL) {
    throw invalid_argument("Could not open image!");
  }
  
  TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &depth);
  depth = (2 << (depth - 1)) - 1;
  
/*  uint16 spp;
  if ((TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp) == 0) || (spp != 1)) {
    throw invalid_argument("Either undefined or unsupported number of samples per pixel");
  }*/
    
  tsize_t stripSize = TIFFStripSize(image);
  int stripMax = TIFFNumberOfStrips(image);
  unsigned long imageOffset = 0;
  
  unsigned long bufferSize = stripMax * stripSize;
  try {
    buffer = new unsigned char[bufferSize];
  } catch (bad_alloc&) {
    throw runtime_error("Failure to allocate storage for tiff image");
  }
  
  for (int stripCount = 0; stripCount < stripMax; ++stripCount) {
    long result = TIFFReadEncodedStrip(image, stripCount, buffer + imageOffset, stripSize);
	if (result == -1) {
	  throw invalid_argument("Read error for tiff image");
	}
    imageOffset += result;
  }
  
  // find the number of rows and number of columns
  if (TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &cols) == 0) {
    throw invalid_argument("Tiff image does not define its width");
  }
  if (TIFFGetField(image, TIFFTAG_IMAGELENGTH, &rows) == 0) {
    throw invalid_argument("Tiff image does not define its length");
  }

  TIFFClose(image);
}
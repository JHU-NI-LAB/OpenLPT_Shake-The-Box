/* 
 *  Tiff.h
 *
 *  A simple wrapper class for reading in tiff files. Tiff uses the libtiff libraries, 
 *  and must be linked against them with -ltiff. 
 *  Tiff inherits from Image.
 *
 *  Last update: 10/2/09 by NTO
 *
 */
 
#ifndef TIFF_H
#define TIFF_H

#include <string>
#include <stdexcept>
#include <iostream>

#include <Image.h>

class Tiffload : public Image {
public:
  // constructor: takes a filename for a tiff image
  Tiffload(std::string f) throw(std::invalid_argument, std::runtime_error);
  // destructor
//  ~Tiffload() {  //buffer will be clear in public object Image
//	  delete[] buffer;
//  };
	
	void Open() throw(std::invalid_argument, std::runtime_error);
		
private:
	std::string filename;
};

#endif // TIFF_H

/*
 *  Image.h
 *  
 *  Base class for various image formats. All image format types should inherit from Image.
 *
 *  Last Update: 10/1/09 by NTO
 *
 */

#ifndef IMAGE_H
#define IMAGE_H

#include <string>
#include <stdexcept>
#include <cstdlib>
#include <iostream>

/*
 * Modified by Shiyong Tan, 1/30/2018
 * memcpy() needs head file: cstring
 * Start:
 */
#include <cstring>
// End

class Image {
public:
	Image();
	Image(int r, int c, int d, unsigned char* pix) throw(std::runtime_error);
	Image(const Image& i) throw(std::runtime_error);
	virtual ~Image();
	
	virtual void Open() throw(std::invalid_argument, std::runtime_error);
	
	int GetRows() const;
	int GetCols() const;
	int GetDepth() const;
	
	// note: pixels_orig must be *pre-allocated*
	void GetPixels(int**& pixels) const throw(std::runtime_error);
	
protected:
	int rows;
	int cols;
	int depth;
	unsigned char* buffer;
};

inline Image::Image() : buffer(NULL)
{}

inline Image::Image(int r, int c, int d, unsigned char* pix) throw(std::runtime_error)
: rows(r), cols(c), depth(d)
{
	if (pix == NULL) {
		throw std::runtime_error("NULL input in Image constructor");
	}
	try {
		buffer = new unsigned char[rows * cols];
	} catch (std::bad_alloc&) {
		throw std::runtime_error("Failure to allocate memory in Image copy constructor");
	}
	memcpy(buffer, pix, rows * cols);
}

inline Image::Image(const Image& i) throw(std::runtime_error)
: rows(i.rows), cols(i.cols), depth(i.depth)
{
	if (i.buffer == NULL) {
		throw std::runtime_error("NULL input in Image copy constructor");
	}
	try {
		buffer = new unsigned char[rows * cols];
	} catch (std::bad_alloc&) {
		throw std::runtime_error("Failure to allocate memory in Image copy constructor");
	}
	memcpy(buffer, i.buffer, rows * cols);
}

inline Image::~Image() 
{
	if (buffer != NULL) {
	  delete []buffer;
	}
}

inline void Image::Open() throw(std::invalid_argument, std::runtime_error)
{
  // nothing to do for generic Image objects, but we don't want an abstract class
}

inline int Image::GetRows() const
{
	return rows;
}

inline int Image::GetCols() const
{
	return cols;
}

inline int Image::GetDepth() const
{
	return depth;
}

#endif // IMAGE_H

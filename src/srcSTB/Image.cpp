/*
 *  Image.cpp
 *  
 *  Implementation file for Image objects.
 *
 *  Last update: 10/1/09 by NTO
 *
 */

#include <Image.h>

using namespace std;

void Image::GetPixels(int**& pixels) const throw(runtime_error)
{
	if (buffer == NULL) {
		throw runtime_error("No Image stored in Image::GetPixels()");
	}
	if (pixels == NULL) {
		throw runtime_error("NULL argument in Image::GetPixels()");
	}
	// remove after testing with matlab image
	for (int j = 0; j < cols ; ++j)
		pixels[0][j] = 0;

	for (int j = 0; j < rows; ++j)
		pixels[j][0] = 0;


//	 edit after matlab testing
	for (int i = 0; i < rows-1; ++i) {
		for (int j = 0; j < cols-1; ++j) {
			pixels[i+1][j+1] = static_cast<int>(buffer[cols * i + j]);
		}
	}

//	for (int i = 0; i < rows; ++i) {
//		for (int j = 0; j < cols; ++j)
//			pixels[i][j] = static_cast<int>(buffer[cols * i + j]);
//	}
}

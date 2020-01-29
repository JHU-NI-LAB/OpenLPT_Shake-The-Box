/*
 *  Camera.h
 *  
 *
 *  Created by Nicholas T. Ouellette on 10/28/11.
 *  Copyright 2011 Yale University. All rights reserved.
 *
 
 *  Revised by Rui Ni on 2/21/2017
 *  add Noffw and Noffh, offset for the center of image, from the calibration file
 *  
 *
 */



#ifndef CAMERA_H
#define CAMERA_H

#include <fstream>

#include <Position.h>
#include <Matrix.h>
#include <string.h>
#include <limits.h>

class Camera {
public:
	Camera(std::istream& is);
	Camera(const Camera& c);
	~Camera() {};
	
	// return camera projective center (in world coordinates)
	Position Center() const;
	
	// remove distortion; return centered coordinates in physical units
	Position UnDistort(const Position& p) const;
	// add distortion back; return normal images coordinates in pixel units
	Position Distort(const Position& p) const;
	
	// project a position (distorted, in pixels_orig!) on the image plane to 3D world coordinates (in mm)
	Position ImageToWorld(const Position& p) const;
	// project a 3D world position (in mm) to a position on the image plane (undistorted, in mm)
	Position WorldToImage(const Position& p) const;

	// Getting the parameters of the camera images
	int Get_Npixh();
	int Get_Npixw();
	double Get_kr();
	double Get_hpix();
	double Get_wpix();
	double Get_Noffh();
	double Get_Noffw();

	//Setting the number of pixels_orig (Just for testing purposes, remove after testing)
	void Set_Npixh(int i);
	void Set_Npixw(int i);


private:
	// model parameters; names should be more or less the same as from calibTsai.m!
    double Noffh;
    double Noffw;
	int Npixw;
	int Npixh;
	double wpix;
	double hpix;
	double f_eff;
	double kr;
	double kx;
	Matrix R;
	Position T;
	Matrix Rinv;
	Position Tinv;
		
};

inline Camera::Camera(const Camera& c)
: Noffh(c.Noffh), Noffw(c.Noffw), Npixw(c.Npixw), Npixh(c.Npixh), wpix(c.wpix), hpix(c.hpix), f_eff(c.f_eff),
kr(c.kr), kx(c.kx), R(c.R), T(c.T), Rinv(c.Rinv), Tinv(c.Tinv)
{}

inline Position Camera::Center() const
{
	return Tinv;
}

inline Position Camera::UnDistort(const Position& p) const
{
	//Position centered(p);
	// shift origin to the center of the image
	//centered -= Position(Npixw/2, Npixh/2, 0);
	
	// account for left-handed coordinate system...
	Position centered(p.X() - Npixw/2 - Noffw, -p.Y() + Npixh/2 - Noffh, p.Z());
	
	// scale into physical units
	centered *= Position(wpix, hpix, 0);
	// remove possible cylindrical distortion
	// 	cout << corrframes[0] << endl;should this be 1.0/kx? i.e., should kx be bigger or smaller than 1?
	centered *= Position(kx, 1, 0);
	// compute the radial correction factor
	
	 // Old code
	
	if (kr == 0) 
		return centered;
	
	else {
		// New Code (accounting for radial distortion)
		double y = centered.Y();
		double x = centered.X();
		double a = centered.X() / centered.Y();

		y = (1 - sqrt(1 - 4 * pow(y, 2) * kr*(pow(a, 2) + 1))) / (2 * y * kr * (pow(a, 2) + 1));
		x = a*y;

		Position M(x, y, 0);
		return M;
	}
}

inline Position Camera::Distort(const Position& p) const
{
	// compute the radial distortion factor

	double rad = 1 + kr * p.Magnitude2();
	Position pixelcoords(p / rad);
	// remove potential cylindrical distortion
	pixelcoords *= Position(1.0 / kx, 1, 0);
	// scale into pixel units
	pixelcoords *= Position(1.0 / wpix, 1.0 / hpix, 0);
	// shift origin
	return Position(pixelcoords.X() + Npixw/2 + Noffw, -1.0 * (pixelcoords.Y() - Npixh/2)-Noffh, 0);
}

inline Position Camera::ImageToWorld(const Position& p) const
{
	Position tmp(p * T.Z() / f_eff);
	Position proj(tmp.X(), tmp.Y(), T.Z());
	// return coordinates of p in 3D world coordinates
	//return ((proj - T) * R);
	return (Rinv * (proj - T));
}

inline Position Camera::WorldToImage(const Position& p) const
{
	//Position proj(p * Rinv + T);
	Position proj(R * p + T);
	return (proj * (f_eff / proj.Z())); 
}

inline int Camera::Get_Npixh()
{
	return Npixh;
}

inline int Camera::Get_Npixw()
{
	return Npixw;
}

inline double Camera::Get_kr()
{
	return kr;
}

inline double Camera::Get_hpix()
{
	return hpix;
}

inline double Camera::Get_wpix()
{
	return wpix;
}

inline double Camera::Get_Noffw()
{
	return Noffw;
}

inline double Camera::Get_Noffh()
{
	return Noffh;
}

inline void Camera::Set_Npixh(int i)
{
	Camera::Npixh = i;
}

inline void Camera::Set_Npixw(int i)
{
	Camera::Npixw = i;
}


#endif // CAMERA_H

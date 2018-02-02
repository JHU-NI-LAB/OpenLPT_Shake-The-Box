/*
 *  Position.h
 *
 *  A Position object holds the coordinates of a single particle, and 
 *  has lots of useful operators defined.
 *
 *  Last update: 10/28/11 by NTO (made 3D)
 *
 */

#ifndef POSITION_H
#define POSITION_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

class Position {
public:
  // constructor: does nothing
  Position() {};
  // constructor: initializes the coordinates
  Position(double NewX, double NewY, double NewZ);
  Position(double NewX, double NewY, double NewZ, double Newx1, double Newy1, double Newx2, double Newy2, double Newx3, double Newy3, double Newx4, double Newy4, double Newinfo);
  // copy-constructor
  Position(const Position& p);
  // destructor: nothing to do
  ~Position() {};
  
  // set this position to be a part of a track
  void SetTracked();
  // check if this position is a part of a track
  bool IsTracked() const;

  // set this position as fake
  void SetFake();
  // check if this position is a fake
  bool IsFake() const;

  // get x component
  double X() const;
  // get y component
  double Y() const;
	// get z component
  double Z() const;
  
  // set x component
  void Set_X(double a);
  // set y component
  void Set_Y(double a);
  // set z component
  void Set_Z(double a);
  // set info
  void Set_info(double a);
  // if the belongs to track, set the type of track
  void Set_trackType(int a);

  double X1() const;
  double Y1() const;
  double X2() const;
  double Y2() const;
  double X3() const;
  double Y3() const;
  double X4() const;
  double Y4() const;
  double Info() const;
  
  // get the squared Euclidean distance between two Positions
  // (don't take the square root for efficiency)
  friend double Distance(const Position& p1, const Position& p2);
  // get the magnitude of the vector
  double Magnitude() const;
  // get the squared magnitude
  double Magnitude2() const;

	// scalar product
  friend double Dot(const Position& left, const Position& right);
	// element-wise multiplication
	friend const Position Multiply(const Position& left, const Position& right);
	
  // member operators
  Position& operator=(const Position& p);
  // vector sum and difference
  Position& operator+=(const Position& right);
  Position& operator-=(const Position& right);
  // scalar multiplication and division
  Position& operator*=(double right);
  Position& operator/=(double right);
  // element by element multiplication
  Position& operator*=(const Position& right);

  // non-member operators

  // comparison
  friend int operator==(const Position& left, const Position& right);
  friend int operator!=(const Position& left, const Position& right);
	// compare y values (for sorting)
	friend int operator<(const Position& left, const Position& right);
	friend int operator>(const Position& left, const Position& right);
  // vector sum and difference
  friend const Position operator+(const Position& left, const Position& right);
  friend const Position operator-(const Position& left, const Position& right);
  // scalar multiplication and division
  friend const Position operator*(const Position& left, double right);
  friend const Position operator*(double left, const Position& right);
  friend const Position operator/(const Position& left, double right);
  // printing
  friend std::ostream& operator<<(std::ostream& os, const Position& p);

	
	
private:
  double x;
  double y;
  double z;
  double x1;
  double y1;
  double x2;
  double y2;
  double x3;
  double y3;
  double x4;
  double y4;
  double info;
  // a flag specifying whether this position is a part of a track already
  bool tracked;  
  // set track type (0->Inactive/Not a part of track, 1->ShortActiveTrack, 2->LongActvieTrack, 3->ExitTrack, 4->LongActiveNewTrack)
  short int trackType;

  // a flag specifying whether this position is a fake
  bool fake;
};

//  Inline Function Definitions

inline Position::Position(double NewX, double NewY, double NewZ) 
: x(NewX), y(NewY), z(NewZ), tracked(false), fake(false)
{}

inline Position::Position(double NewX, double NewY, double NewZ, double Newx1, double Newy1, double Newx2, double Newy2, double Newx3, double Newy3, double Newx4, double Newy4, double Newinfo) 
: x(NewX), y(NewY), z(NewZ), x1(Newx1), y1(Newy1), x2(Newx2), y2(Newy2), x3(Newx3), y3(Newy3), x4(Newx4), y4(Newy4), info(Newinfo), tracked(false), fake(false)
{}

inline Position::Position(const Position& p) 
: x(p.x), y(p.y), z(p.z), x1(p.x1), y1(p.y1), x2(p.x2), y2(p.y2), x3(p.x3), y3(p.y3), x4(p.x4), y4(p.y4), info(p.info), tracked(p.tracked), fake(p.fake)
{}

inline void Position::SetTracked()
{
  tracked = true;
}

inline bool Position::IsTracked() const
{
  return tracked;
}

inline void Position::SetFake()
{
	fake = true;
}

inline bool Position::IsFake() const
{
	return fake;
}

inline double Position::Magnitude() const
{
  return sqrt(x * x + y * y + z * z);
}

inline double Position::Magnitude2() const
{
  return (x * x + y * y + z * z);
}

inline double Position::X() const
{
	return x;
}

inline double Position::Y() const
{
	return y;
}

inline double Position::Z() const
{
	return z;
}

inline void Position::Set_X(double a) 
{
	Position::x = a;
}

inline void Position::Set_Y(double a) 
{
	Position::y = a;
}

inline void Position::Set_Z(double a)
{
	Position::z = a;
}
inline void Position::Set_info(double a)
{
	Position::info = a;
}

inline void Position::Set_trackType(int a)
{
	trackType = a;
}

inline double Position::X1() const
{
	return x1;
}

inline double Position::Y1() const
{
	return y1;
}


inline double Position::X2() const
{
	return x2;
}

inline double Position::Y2() const
{
	return y2;
}


inline double Position::X3() const
{
	return x3;
}

inline double Position::Y3() const
{
	return y3;
}


inline double Position::X4() const
{
	return x4;
}

inline double Position::Y4() const
{
	return y4;
}

inline double Position::Info() const
{
	return info;
}

inline Position& Position::operator=(const Position& p)
{
  x = p.x;
  y = p.y;
  z = p.z;
  x1 = p.x1;
  y1 = p.y1;
  x2 = p.x2;
  y2 = p.y2;
  x3 = p.x3;
  y3 = p.y3;
  x4 = p.x4;
  y4 = p.y4;
  info = p.info;
  tracked = p.tracked;
  fake = p.fake;

  return *this;
}

inline Position& Position::operator+=(const Position& right)
{
  x += right.x;
  y += right.y;
	z += right.z;

  return *this;
}

inline Position& Position::operator-=(const Position& right)
{
  x -= right.x;
  y -= right.y;
	z -= right.z;

  return *this;
}

inline Position& Position::operator*=(double right)
{
  x *= right;
  y *= right;
	z *= right;

  return *this;
}

inline Position& Position::operator/=(double right)
{
  x /= right;
  y /= right;
	z /= right;

  return *this;
}

inline Position& Position::operator*=(const Position& right)
{
	x *= right.x;
	y *= right.y;
	z *= right.z;
	
	return *this;
}


#endif // POSITION_H

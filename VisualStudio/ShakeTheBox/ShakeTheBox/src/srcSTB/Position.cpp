/*
 *  Position.cpp
 *
 *  Implementation file for Position objects.
 *
 *  Last update: 10/28/11 by NTO (made 3D)
 *
 */

#include <Position.h>

using namespace std;

double Distance(const Position& p1, const Position& p2) 
{
  return ((p1.x - p2.x) * (p1.x - p2.x) 
					+ (p1.y - p2.y) * (p1.y - p2.y)
					+ (p1.z - p2.z) * (p1.z - p2.z));
}

int operator==(const Position& left, const Position& right) 
{
  return ((left.x == right.x) 
					&& (left.y == right.y) 
					&& (left.z == right.z));
}

int operator!=(const Position& left, const Position& right) 
{
  return ((left.x != right.x) 
					|| (left.y != right.y)
					|| (left.z != right.z));
}

int operator<(const Position& left, const Position& right)
{
	return left.y < right.y;
}

int operator>(const Position& left, const Position& right)
{
	return left.y > right.y;
}

const Position operator+(const Position& left, const Position& right)
{
  return Position(left.x + right.x, left.y + right.y, left.z + right.z);
}

const Position operator-(const Position& left, const Position& right)
{
  return Position(left.x - right.x, left.y - right.y, left.z - right.z);
}

const Position operator*(const Position& left, double right)
{
  return Position(left.x * right, left.y * right, left.z * right);
}

const Position operator*(double left, const Position& right)
{
  return Position(left * right.x, left * right.y, left * right.z);
}

const Position operator/(const Position& left, double right)
{
  return Position(left.x / right, left.y / right, left.z / right);
}

double Dot(const Position& left, const Position& right) 
{
  return ((left.x * right.x) + (left.y * right.y) + (left.z * right.z));
}

const Position Multiply(const Position& left, const Position& right)
{
	return Position(left.x * right.x, left.y * right.y, left.z * right.z);
}

ostream& operator<<(ostream& os, const Position& p)
{
	os << "[" << p.x << " " << p.y << " " << p.z << "] [" << p.x1 << " " << p.y1 << "] [" << p.x2 << " " << p.y2 << "] [" << p.x3 << " " << p.y3 << "] [" << p.x4 << " " << p.y4 << "] [" << p.info << "\t";
	if (p.fake) {
		os << 1 << "]";
	} else {
		os << 0 << "]";
	}

  return os;
}

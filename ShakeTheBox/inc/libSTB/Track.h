/*
 * ----------------------------------------------------------------------------
 * File: Track.h
 * ----------------------------------------------------------------------------
 * This is the header file for Track objects.  A Track object is essentially
 * a fancy frame: it inherits from Frame, but contains more information.
 * ----------------------------------------------------------------------------
 * Created 7/17/03
 * Last updated 3/2/04
 * ----------------------------------------------------------------------------
 * Modified for Yale by NTO on 10/28/11
 *	No longer inherits from Frame
 * ----------------------------------------------------------------------------
 */

#ifndef TRACK_H
#define TRACK_H

#include <deque>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <assert.h>

#include <Position.h>

class Track {

public:

  // default constructor
  Track();
  // constructor: give one Position argument, and a time (defaults to 0)
  Track(const Position& p, int t = 0);
  // copy-constructor
  Track(const Track& t);
  // destructor: nothing to do
  ~Track() {};

  // Add a new point to the Track
  void AddNext(const Position& p, int t);
  // Add another Track onto the end of this one
  void AddNext(const Track& t);

  // Delete a point to the end of the Track
  void DeleteBack();

  void DeleteFront();

  // Add a new point in the begining Track
  void AddFront(const Position& p, int t);
  // Add another Track onto the begining of this one
  void AddFront(const Track& t);

  const Position First() const;
  const Position Second() const;
  const Position Third() const;

  // return the last point on the Track
  const Position Last() const;
  // return the penultimate point on the Track
  const Position Penultimate() const;
  // return the antepenultimate point on the Track
  const Position Antepenultimate() const;
  

  // get the length of the Track (not including any ending extrapolated points)
  int Length() const;
  // get the time of a particular element
  int GetTime(int index) const throw(std::out_of_range);
  // get the position of nth element
  const Position GetPos(int n) const throw(std::out_of_range);;

  // check the occlusion counter
  int OcclusionCount() const;
  // increment the occlusion counter
  void Occluded();
  // reset the occlusion counter
  void ResetCounter();

  // find the number of usable fake points on the track
  int NumFake() const;

  // write the track as part of a GDF file (see Tracker.h for format)
  void WriteGDF(std::ofstream& output, float index, float fps = 1) const;
  
  // member operators
  Track& operator=(const Track& t);
  // return the nth point on the Track
  const Position operator[](int n);
  bool operator ==(Track& t1);

  // non-member operators
  friend std::ostream& operator<<(std::ostream& os, const Track& t);

  // print only the estimated points
  void PrintEstimates(std::ostream& os) const;
  
  // set track as active / inactive
  void SetActive(bool activity) { (active = activity); }
  
  // check if the track is active / inactive
  bool IsActive() { return active; }

  // check if the track has a particle at time 't'
  bool Exists(int t);

  // free memory
  void Clear();

private:
	std::deque<Position> pos;
  std::deque<int> time;	// the time (as an integer frame number)
  int occluded; // a counter keeping track of the number of frames this 
                // track hasn't had a real particle added to it.
  int npoints;

  bool active;
};

// Inline Function Definitions

inline Track::Track() : occluded(0), npoints(0) { active = true; }

inline Track::Track(const Position& p, int t /* = 0 */) 
: occluded(0), npoints(1)
{
	pos.push_back(p);
  time.push_back(t);
  assert(active = true);
}

inline Track::Track(const Track& t)
: pos(t.pos), time(t.time), occluded(t.occluded), npoints(t.npoints)
{}

inline void Track::AddNext(const Position& p, int t)
{
  pos.push_back(p);
  time.push_back(t);
  assert(active = true);
  ++npoints;
}

inline void Track::DeleteBack() {
	pos.pop_back();
	time.pop_back();
	--npoints;
}

inline void Track::DeleteFront() {
	pos.pop_front();
	time.pop_front();
	--npoints;
}

inline void Track::AddFront(const Position& p, int t)
{
	pos.push_front(p);
	time.push_front(t);
	//assert(active = true);
	++npoints;
}

inline const Position Track::First() const
{
  return pos[0];
}

inline const Position Track::Second() const
{
  return pos[1];
}

inline const Position Track::Third() const
{
  return pos[2];
}

inline const Position Track::Last() const
{
  return pos[npoints - 1];
}

inline const Position Track::Penultimate() const
{
  return pos[npoints - 2];
}

inline const Position Track::Antepenultimate() const
{
  return pos[npoints - 3];
}

inline int Track::GetTime(int index) const throw(std::out_of_range)
{
  try {
    return time.at(index);
  } catch (std::out_of_range& e) {
    std::cerr << e.what() << std::endl;
    throw std::out_of_range("Caught out_of_range in Track::GetTime()");
  }
}

inline const Position Track::GetPos(int n) const throw(std::out_of_range)
{
	try {
		return pos.at(n);
	}
	catch (std::out_of_range& e) {
		std::cerr << e.what() << std::endl;
		throw std::out_of_range("Caught out_of_range in Track::GetPos()");
	}
}

inline int Track::OcclusionCount() const
{
  return occluded;
}

inline void Track::Occluded() 
{
  ++occluded;
}

inline void Track::ResetCounter()
{
  occluded = 0;
}

inline Track& Track::operator=(const Track& t)
{
  pos = t.pos;
  time = t.time;
  occluded = t.occluded;
  npoints = t.npoints;

  return *this;
}

inline const Position Track::operator[](int n)
{
	return pos[n];
}

inline bool Track::operator ==(Track& t1) {
	if (t1.GetPos(0) == pos[0] && t1.Last() == pos[npoints - 1])
		return true;
}

inline void Track::Clear()
{
  pos.clear();
  time.clear();
  npoints = 0;
  occluded = 0;
}

inline bool Track::Exists(int t) {
	bool exists = false;
	if (t <= time[npoints - 1] && t >= time[0])
		exists = true;

	return exists;
}
#endif /* TRACK_H */

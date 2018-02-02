/*
 *  Frame.h
 *
 *  A Frame holds the locations of all the located particles at a given time step.
 *
 *  Last update: 10/6/09 by NTO
 *
 */

#ifndef FRAME_H
#define FRAME_H

#include <string>
#include <deque>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <Position.h>

class Frame {
public:
	// simple constructor: do nothing, and wait for the user to give us data
	Frame();
	// constructor: generate a Frame with one point
	Frame(const Position& p);
	// constructor: create a Frame from a deque of Positions
	Frame(const std::deque<Position>& p);
	// copy-constructor
	Frame(const Frame& f);
	// destructor: nothing to do
	//~Frame();

	// get the number of particles in the Frame
	int NumParticles() const;

	// get a certain position from the deque
	Position Get_pos(int i);

	// free the memory associated with this Frame
	void Clear();

	// fill an already-constructed Frame with data
	void Fill(const std::deque<Position>& p);

	// add an extra position to the deque
	void Add(const Position& p);

	// delete a Position at index i 
	void Delete(int i);

	// member operators

	Frame& operator=(const Frame& f);
	// random access to positions in the Frame
	Position& operator[](const int x) throw(std::out_of_range);

	// non-member operators

	friend std::ostream& operator<<(std::ostream& os, const Frame& f);

	// a const_iterator class: provided for use with the algorithm header
	class const_iterator;
	friend class const_iterator;
	class const_iterator : public std::iterator<std::random_access_iterator_tag, Position, std::ptrdiff_t> {
	public:
		const_iterator();
		const_iterator(const Frame& f);
		const_iterator(const Frame& f, bool);
		const_iterator(const const_iterator& i);

		int where() const;

		bool operator==(const const_iterator& i) const;
		bool operator!=(const const_iterator& i) const;
		bool operator<(const const_iterator& i) const;
		bool operator>(const const_iterator& i) const;
		bool operator<=(const const_iterator& i) const;
		bool operator>=(const const_iterator& i) const;

		const_iterator& operator++();
		const_iterator operator++(int);
		const_iterator& operator--();
		const_iterator operator--(int);


		const_iterator operator+(int n) const;
		const_iterator operator-(int n) const;
		const_iterator& operator+=(int n);
		const_iterator& operator-=(int n);

		int operator-(const const_iterator& i) const;

		const Position& operator*() const;
		const Position* operator->() const;
		const Position& operator[](int n) const;

	private:
		std::deque<Position>::const_iterator it;
		int index;
	};

	const_iterator begin() const;
	const_iterator end() const;

	std::deque<Position>& Get_PosDeque();

protected:
	// the positions of the particles in this frame
	std::deque<Position> pos;

private:
	int nparticles;

};

inline Frame::Frame() : nparticles(0)
{}

inline Frame::Frame(const Position& p)
{
	pos.push_back(p);
	nparticles = 1;
}

inline Frame::Frame(const std::deque<Position>& p) : pos(p)
{
	nparticles = pos.size();

}

inline Frame::Frame(const Frame& f) : pos(f.pos), nparticles(f.nparticles)
{}
	
inline int Frame::NumParticles() const
{
  return nparticles;
}

inline void Frame::Clear() 
{
  pos.clear();
	nparticles = 0;
}

inline void Frame::Fill(const std::deque<Position>& p)
{
	pos = p;
	nparticles = p.size();
}

inline void Frame::Add(const Position& p)
{
	pos.push_back(p);
	nparticles = pos.size();
}

inline void Frame::Delete(int i) {
	pos.erase(pos.begin() + i);
	nparticles = pos.size();

}

inline Frame& Frame::operator=(const Frame& f)
{
  pos = f.pos;
	nparticles = f.nparticles;

  return *this;
}

inline Position& Frame::operator[](const int x) throw(std::out_of_range)
{
	try {
    return pos.at(x);
	} catch (std::out_of_range& e) {
		std::cerr << e.what() << std::endl;
		throw std::out_of_range("Caught out_of_range in Frame::operator[]");
	}
}

inline std::deque<Position>& Frame::Get_PosDeque() {
	return pos;
}

inline Frame::const_iterator::const_iterator() : index(-1)
{
}

inline Frame::const_iterator::const_iterator(const Frame& f) : index(0)
{
	it = f.pos.begin();
}

inline Frame::const_iterator::const_iterator(const Frame& f, bool) : index(-1)
{
	it = f.pos.end();
}

inline Frame::const_iterator::const_iterator(const const_iterator& i) : it(i.it), index(i.index)
{}

inline int Frame::const_iterator::where() const
{
	return index;
}

inline bool Frame::const_iterator::operator==(const const_iterator& i) const
{
  return it == i.it;
}

inline bool Frame::const_iterator::operator!=(const const_iterator& i) const
{
	return !(*this == i);
}

inline bool Frame::const_iterator::operator<(const const_iterator& i) const
{
	return (it < i.it);
}

inline bool Frame::const_iterator::operator>(const const_iterator& i) const
{
	return (it > i.it);
}

inline bool Frame::const_iterator::operator<=(const const_iterator& i) const
{
	return (it <= i.it);
}

inline bool Frame::const_iterator::operator>=(const const_iterator& i) const
{
	return (it >= i.it);
}

inline Frame::const_iterator& Frame::const_iterator::operator++()
{
	++it;
	++index;
	return *this;
}

inline Frame::const_iterator Frame::const_iterator::operator++(int)
{
	const_iterator tmp = *this;
	++*this;
	++index;
	return tmp;
}

inline Frame::const_iterator& Frame::const_iterator::operator--()
{
	--it;
	--index;
	return *this;
}

inline Frame::const_iterator Frame::const_iterator::operator--(int)
{
	const_iterator tmp = *this;
	--*this;
	--index;
	return tmp;
}

inline Frame::const_iterator Frame::const_iterator::operator+(int n) const
{
	const_iterator tmp = *this;
	return tmp += n;
}

inline Frame::const_iterator Frame::const_iterator::operator-(int n) const
{
	const_iterator tmp = *this;
	return tmp -= n;
}

inline Frame::const_iterator& Frame::const_iterator::operator+=(int n) 
{
	it += n;
	index += n;
	return *this;
}

inline Frame::const_iterator& Frame::const_iterator::operator-=(int n)
{
	return *this += -n;
}

inline int Frame::const_iterator::operator-(const const_iterator& i) const
{
	return (it - i.it);
}

inline const Position& Frame::const_iterator::operator*() const
{
	return *it;
}

inline const Position* Frame::const_iterator::operator->() const
{
	return &(*it);
}

inline const Position& Frame::const_iterator::operator[](int n) const
{
	return *(*this + n);
}

inline Frame::const_iterator Frame::begin() const
{
	return const_iterator(*this);
}

inline Frame::const_iterator Frame::end() const
{
	return const_iterator(*this, true);
}


#endif // FRAME_H 

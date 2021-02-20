/*
 *  Frame.cpp
 *
 *  Implementation for Frame objects.
 *
 *  Last update: 10/1/09 by NTO
 *
 */

#include <fstream>
#include <math.h>
#include <Frame.h>


using namespace std;

ostream& operator<<(ostream& os, const Frame& f)
{
  for (unsigned int i = 0; i < f.pos.size(); ++i) {
    os << "\t" << f.pos[i] << "\n";
  }

  return os;
}



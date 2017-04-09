
/* data structures to describe structures in Cartesian space (like points, boxes, etc. */

// COMPILATION AND RUN COMMANDS:
// -
// -

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "basic.h"
#include <vector>
#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<int DIM> struct Point {
  /* represent a point in DIM dimensions */
    
  double x[DIM];              // coordinates
  
  Point (const Point &p) {     // copy constructor
      for(int i = 0; i < DIM; i++) x[i] = p.x[i];
  }
  
  Point& operator= (const Point &p){     // assignment operator
      for(int i = 0; i < DIM; i++) x[i] = p.x[i];
      return *this;
  }
  
  bool operator== (const Point &p) const {
      for(int i = 0; i < DIM; i++) if (x[i] != p.x[i]) return false;
      return true;
  }
  
  Point (double x0 = 0.0, double x1 = 0.0, double x2 = 0.0) {   // constructor by coordinate values
      x[0] = x0;
      if (DIM > 1) x[1] = x1;
      if (DIM > 2) x[2] = x2;
      if (DIM > 3) throw("Point not implemented for DIM > 3");
  }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<int DIM> double dist (const Point<DIM> &p, const Point<DIM> &q) {
  /* return the distance between two points in DIM dimensions */
  
  double dd = 0.0;
  for(int j = 0; j < DIM; j++) dd += (q.x[j]-p.x[j])*(q.x[j]-p.x[j]);
  return sqrt(dd);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<int DIM> double dist_xy (const Point<DIM> &p, const Point<DIM> &q, int i) {
  /* return the distance between two pints in DIM dimensions */
  
  return q.x[i]-p.x[i];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template<int DIM> struct Box {
  /* represent a Cartesian box in DIM dimension */
  
  Point<DIM> lo, hi;      // Diagonally opposite corners (min of all coordinates and max of all coordinates) are stored as two points
  Box() {}    // constructor
  Box(const Point<DIM> &mylo, const Point<DIM> &myhi) : lo(mylo), hi(myhi) {}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// If point p lies outside box b, the distance to the nearest point on b is returned. If p is inside b or on its surface, zero is returned.
template<int DIM> double dist (const Box<DIM> &b, const Point<DIM> &p) {
  /* if point p lies outside box b, the distance to the nearest point on b is returned
  If p is inside b or on its surface, zero is returned
  */
  
  double dd = 0.;
  for(int i = 0; i < DIM; i++) {
      if (p.x[i]<b.lo.x[i]) dd += (p.x[i]-b.lo.x[i])*(p.x[i]-b.lo.x[i]);
      if (p.x[i]>b.hi.x[i]) dd += (p.x[i]-b.hi.x[i])*(p.x[i]-b.hi.x[i]);
  }
  return sqrt(dd);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Node in a binary tree of boxes containing points.
template<int DIM> struct Boxnode : Box<DIM> {
  /* node in a binary tree of boxes containing points */
  
  int mom, dau1, dau2, ptlo, pthi;
  Boxnode() {}
  Boxnode(Point<DIM> mylo, Point<DIM> myhi, int mymom, int myd1, int myd2, int myptlo, int mypthi) : Box<DIM>(mylo, myhi), mom(mymom), dau1(myd1), dau2(myd2), ptlo(myptlo), pthi(mypthi) {}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////

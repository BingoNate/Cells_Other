
// +++++++
// ++ Adapted from Numerical Recipes 3
// +++++++

#pragma once

#include "pointbox.h"
#include "basic.h"
#include "hash.h"
#include <vector>
#include <cmath>
#include <set>
#include <iostream>

using namespace std;


typedef unsigned long long int Ullong;    // 64 bit integer


// Structure of a circle
struct Circle {
    Point<2> center;
    double radius;
    Circle(const Point<2> &cen, double rad) : center(cen), radius(rad) {}
};

Circle circumcircle(Point<2> a, Point<2> b, Point<2> c) {
    double a0,a1,c0,c1,det,asq,csq,ctr0,ctr1,rad2;
    a0 = a.x[0] - b.x[0]; a1 = a.x[1] - b.x[1];
    c0 = c.x[0] - b.x[0]; c1 = c.x[1] - b.x[1];
    det = a0*c1 - c0*a1;
    if (det == 0.0) throw("no circle through colinear points");
    det = 0.5/det;
    asq = a0*a0 + a1*a1;
    csq = c0*c0 + c1*c1;
    ctr0 = det*(asq*c1 - csq*a1);
    ctr1 = det*(csq*a0 - asq*c0);
    rad2 = ctr0*ctr0 + ctr1*ctr1;
    return Circle(Point<2>(ctr0 + b.x[0], ctr1 + b.x[1]), sqrt(rad2));
}


// Structure of triangle elements
struct Triel {
    Point<2> *pts;              // pointer to the array of points
    int p[3];                   // triangle's vertices in counterclockwise order
    int d[3];                   // pointers for up to 3 daughters
    int stat;                   // nonzero if this element is live
    
    // Set data in the triangle
    void setme(int a, int b, int c, Point<2> *ptss) {
        pts = ptss;
        p[0] = a; p[1] = b; p[2] = c;
        d[0] = d[1] = d[2] = -1;
        stat = 1;
    }
    
    // Return 1 if the point is in the triangle, 0 on the boundary, -1 if it's ouside (CCW is assumed)
    int contains(Point<2> point) {
        double d;
        int i,j,ztest = 0;
        for(i=0; i<3; i++) {
            j = (i+1) % 3;
            d = (pts[p[j]].x[0]-pts[p[i]].x[0])*(point.x[1]-pts[p[i]].x[1]) - (pts[p[j]].x[1]-pts[p[i]].x[1])*(point.x[0]-pts[p[i]].x[0]);
            if (d < 0.0) return -1;
            if (d == 0.0) ztest = 1;
        }
        return (ztest ? 0 : 1);
    }
};


// Return positive, zero, or negative if point d is in, on, or outside the circle through points a, b anc c, respectively
double incircle(Point<2> d, Point<2> a, Point<2> b, Point<2> c) {
    Circle cc = circumcircle(a,b,c);
    double radd = (d.x[0]-cc.center.x[0])*(d.x[0]-cc.center.x[0]) + (d.x[1]-cc.center.x[1])*(d.x[1]-cc.center.x[1]);
    return ((cc.radius)*(cc.radius) - radd);
}



// Null hash function. Use a key (assumed to be already hashed) as its own hash.
struct Nullhash {
    Nullhash(int nn) {}
    inline Ullong fn(const void *key) const { return *((Ullong *)key); }
};


// Structure for constructing Delaunay triangulation
struct Delaunay {
    int npts,ntri,ntree,ntreemax,opt;               // Number of points, triangles, elements in Triel list, maximum of the same
    double delx,dely;
    vector<Point<2> > pts;
    vector<Triel> thelist;                          // the list of Triel elements
    Hash<Ullong,int,Nullhash> *linehash;
    Hash<Ullong,int,Nullhash> *trihash;
    vector<int> perm;
    Delaunay(vector<Point<2> > &pvec, int options = 0);
    Ranhash hashfn;                                 // the raw hash function
    double interpolate(const Point<2> &p, const vector<double> &fnvals, double defaultval=0.0);
    void insertapoint(int r);
    int whichcontainspt(const Point<2> &p, int strict = 0);
    int storetriangle(int a, int b, int c);
    void erasetriangle(int a, int b, int c, int d0, int d1, int d2);
    static unsigned int jran;                       // random number counter
    static const double fuzz, bigscale;
};
const double Delaunay::fuzz  = 1.0e-6;
const double Delaunay::bigscale = 1000.0;
unsigned int Delaunay::jran = 14921620;


// Constructor
Delaunay::Delaunay(vector< Point<2> > &pvec, int options) : npts(pvec.size()), ntri(0), ntree(0), ntreemax(10*npts+1000), opt(options), pts(npts+3), thelist(ntreemax) {
    // If bit 0 in options is nonzero, hash memories used in the construction are deleted.
    
    int j;
    double xl,xh,yl,yh;
    
    linehash = new Hash<Ullong,int,Nullhash>(6*npts+12,6*npts+12);
    trihash = new Hash<Ullong,int,Nullhash>(2*npts+6,2*npts+6);
    
    //perm = new int[npts];       // permutation for randomizing point order
    //perm.reserve(npts);
    
    // Choose the initial triangle large enough to cover all points, but these points are fictitious
    xl = xh = pvec[0].x[0];
    yl = yh = pvec[0].x[1];
    
    // Copy points to local store and calculate their bounding box
    for (j=0; j<npts; j++) {
        pts[j] = pvec[j];
        perm.push_back(j);
        if (pvec[j].x[0] < xl) xl = pvec[j].x[0];
        if (pvec[j].x[0] > xh) xh = pvec[j].x[0];
        if (pvec[j].x[1] < yl) yl = pvec[j].x[1];
        if (pvec[j].x[1] > yh) yh = pvec[j].x[1];
    }
    delx = xh - xl;
    dely = yh - yl;
    
    // Initial points
    pts[npts] = Point<2>(0.5*(xl + xh), yh + bigscale*dely);
    pts[npts+1] = Point<2>(xl - 0.5*bigscale*delx,yl - 0.5*bigscale*dely);
    pts[npts+2] = Point<2>(xh + 0.5*bigscale*delx,yl - 0.5*bigscale*dely);
    storetriangle(npts,npts+1,npts+2);
    
    for (j=npts; j>0; j--) SWAP(perm[j-1],perm[hashfn.int64(jran++) % j]);
    
    // Real deal is going on here
    for (j=0; j<npts; j++) insertapoint(perm[j]);
    
    // Delete the huge root triangle and all its connecting edges
    for (j=0; j<ntree; j++) {
        if (thelist[j].stat > 0) {
            if (thelist[j].p[0] >= npts || thelist[j].p[1] >= npts ||
                thelist[j].p[2] >= npts) {
                thelist[j].stat = -1;
                ntri--;
            }
        }
    }
    if (!(opt & 1)) {
        perm.clear();
        delete trihash;
        delete linehash;
    }
}


// Insert a point to the triangulation
void Delaunay::insertapoint(int r) {
    
    int i,j,k,l,s,tno,ntask,d0,d1,d2;
    Ullong key;
    
    int tasks[50], taski[50], taskj[50];        // stacks
    
    // Find the triangle containing point. Fuzz it if it's on the edge.
    for (j=0; j<3; j++) {
        tno = whichcontainspt(pts[r],1);
        if (tno >= 0) break;                // the desired result, point is okay
        pts[r].x[0] += fuzz * delx * (hashfn.doub(jran++)-0.5);
        pts[r].x[1] += fuzz * dely * (hashfn.doub(jran++)-0.5);
    }
    
    if (j == 3) throw("points degenerate even after fuzzing");
    
    ntask = 0;
    
    // Points in the triangle
    i = thelist[tno].p[0]; j = thelist[tno].p[1]; k = thelist[tno].p[2];
    
    // The following line is used by convex hull application, it causes any points already known to be interior to the convex hull to be omitted from the triangulation
    if (opt & 2 && i < npts && j < npts && k < npts) return;
    
    // Create 3 triangles and queue them in the stacks
    d0 = storetriangle(r,i,j);
    tasks[++ntask] = r; taski[ntask] = i; taskj[ntask] = j;
    d1 = storetriangle(r,j,k);
    tasks[++ntask] = r; taski[ntask] = j; taskj[ntask] = k;
    d2 = storetriangle(r,k,i);
    tasks[++ntask] = r; taski[ntask] = k; taskj[ntask] = i;
    
    // Erase the old triangle
    erasetriangle(i,j,k,d0,d1,d2);
    
    // Legalize edges recursively
    while (ntask) {
        
        s=tasks[ntask]; i=taski[ntask]; j=taskj[ntask--];
        
        key = hashfn.int64(j) - hashfn.int64(i);                // look up fourth point
        
        if ( ! linehash->get(key,l) ) continue;                 // case of no triangle on other side
        
        // Needs legalizing? Do an edge flip.
        if (incircle(pts[l],pts[j],pts[s],pts[i]) > 0.0){
            
            // Create two new triangles
            d0 = storetriangle(s,l,j);
            d1 = storetriangle(s,i,l);
            
            // And erase the old ones
            erasetriangle(s,i,j,d0,d1,-1);
            erasetriangle(l,j,i,d0,d1,-1);
            
            // Erase line in both directions
            key = hashfn.int64(i)-hashfn.int64(j);
            linehash->erase(key);
            key = 0 - key;              // unsigned, hence binary minus
            linehash->erase(key);
            
            // Two new edges now need checking
            tasks[++ntask] = s; taski[ntask] = l; taskj[ntask] = j;
            tasks[++ntask] = s; taski[ntask] = i; taskj[ntask] = l;
        }
    }
}


// Given point p, return index in 'thelist' of the triangle in the triangulation that containts it, or return -1 for failure. If strict is nonzero, require strict containment, otherwise allow the point to lie on the edge.
int Delaunay::whichcontainspt(const Point<2> &p, int strict) {
    
    int i,j,k=0;
    
    while (thelist[k].stat <= 0) {
        
        for (i=0; i<3; i++) {
            
            if ((j = thelist[k].d[i]) < 0) continue;
            if (strict) {
                if (thelist[j].contains(p) > 0) break;
            } else {
                if (thelist[j].contains(p) >= 0) break;
            }
            
        }
        
        if (i == 3) return -1;
        k = j;
    }
    
    return k;
}


// Erase triangle abc in trihash and inactivate it in 'thelist' after setting its daughters
void Delaunay::erasetriangle(int a, int b, int c, int d0, int d1, int d2) {
    
    Ullong key;
    int j;
    key = hashfn.int64(a) ^ hashfn.int64(b) ^ hashfn.int64(c);
    if (trihash->get(key,j) == 0) throw("nonexistent triangle");
    trihash->erase(key);
    thelist[j].d[0] = d0; thelist[j].d[1] = d1; thelist[j].d[2] = d2;
    thelist[j].stat = 0;
    ntri--;
    
}


// Store a triangle with vertices a,b,c in trihash. Store its points in linehash under keys to opposite sited. Add it to 'thelist', returning its index there.
int Delaunay::storetriangle(int a, int b, int c) {
    
    Ullong key;
    thelist[ntree].setme(a,b,c,&pts[0]);
    key = hashfn.int64(a) ^ hashfn.int64(b) ^ hashfn.int64(c);
    trihash->set(key,ntree);
    key = hashfn.int64(b)-hashfn.int64(c);
    linehash->set(key,a);
    key = hashfn.int64(c)-hashfn.int64(a);
    linehash->set(key,b);
    key = hashfn.int64(a)-hashfn.int64(b);
    linehash->set(key,c);
    if (++ntree == ntreemax) throw("thelist is sized too small");
    ntri++;
    return (ntree-1);
    
}


// Triangular grid interpolation of a function. Given an arbitrary point p and a vector of function values fnvals at the points that were used to construct the Delaunay structure, return the linearly interpolated function value in the triangle in which p lies. If p lies outside of the triangulation, instead return defaultval.
double Delaunay::interpolate(const Point<2> &p, const vector<double> &fnvals, double defaultval) {
    int n,i,j,k;
    double wgts[3];
    double ipts[3];
    double sum, ans = 0.0;
    n = whichcontainspt(p);
    if (n < 0) return defaultval;
    for (i=0; i<3; i++) ipts[i] = thelist[n].p[i];
    for (i=0,j=1,k=2; i<3; i++,j++,k++) {
        if (j == 3) j=0;
        if (k == 3) k=0;
        wgts[k]=(pts[ipts[j]].x[0]-pts[ipts[i]].x[0])*(p.x[1]-pts[ipts[i]].x[1]) - (pts[ipts[j]].x[1]-pts[ipts[i]].x[1])*(p.x[0]-pts[ipts[i]].x[0]);
    }
    sum = wgts[0] + wgts[1] + wgts[2];
    if (sum == 0) throw("degenerate triangle");
    for (i=0; i<3; i++) ans += wgts[i]*fnvals[ipts[i]]/sum;
    return ans;
}


// Structure for an edge in a Voronoi diagram, containing its two endpoints and an integer pointer to the site of which it is a boundary
struct Voredge {
    Point<2> p[2];
    int nearpt;
    Voredge() {}
    Voredge(Point<2> pa, Point<2> pb, int np) : nearpt(np) {
        p[0] = pa; p[1] = pb;
    }
};


struct Voronoi : Delaunay {
    int nseg;                           // number of edges in the diagram
    vector<int> trindx;                 // will index triangles
    vector<Voredge> segs;               // will be array of all segments
    Voronoi(vector<Point<2> > pvec);    // constructor from array of points
};


// Constructor
Voronoi::Voronoi(vector<Point<2> > pvec) :
Delaunay(pvec,1), nseg(0), trindx(npts), segs(6*npts+12) {
    int i,j,k,p,jfirst;
    Ullong key;
    Triel tt;
    Point<2> cc, ccp;
    
    // Create a table so that given a point number, we can find one triangle with it as a vertex. Wrap a 3D triangle to a 1D point list.
    for (j=0; j<ntree; j++) {
        if (thelist[j].stat <= 0) continue;
        tt = thelist[j];
        for (k=0; k<3; k++) trindx[tt.p[k]] = j;
    }
    
    // Now loop over the sites
    for (p=0; p<npts; p++) {
        
        tt = thelist[trindx[p]];
        if (tt.p[0] == p) {i = tt.p[1]; j = tt.p[2];}           // Get the vertices into canonical order
        else if (tt.p[1] == p) {i = tt.p[2]; j = tt.p[0];}
        else if (tt.p[2] == p) {i = tt.p[0]; j = tt.p[1];}
        else throw("triangle should contain p");
        
        // Save starting index and its circumcircle
        jfirst = j;
        ccp = circumcircle(pts[p],pts[i],pts[j]).center;
        
        // Go around counterclockwise, find circumcenters and store segments
        while (1) {
            key = hashfn.int64(i) - hashfn.int64(p);
            if ( ! linehash->get(key,k) ) throw("Delaunay is incomplete");
            cc = circumcircle(pts[p],pts[k],pts[i]).center;
            segs[nseg++] = Voredge(ccp,cc,p);
            if (k == jfirst) break;         // circumnavigation completed, normal way out
            ccp = cc;
            j=i;
            i=k;
        }
    }
}


// Vertex
struct Vertex{
    Vertex(int idx) : idx(idx) {}
    int idx;                        // vertex index
    set<int> neighs;                // neighbour list
    void addNeigh( int a );
};

void Vertex::addNeigh( int a ){
    
    neighs.insert( a );
    
}

// Neighbor list
struct Neighbour : Delaunay {
    int vrtx_size;
    vector<Vertex> vertices;
    vector<int> trindx;
    double threshold;                                                                           // threshold value for determination of good edges
    Neighbour(vector<Point<2> > pvec, double thr, int m_size);                                  // constructor
};

Neighbour::Neighbour(vector<Point<2> > pvec, double thr, int m_size) : Delaunay(pvec,1), threshold(thr*thr), trindx(npts), vrtx_size(m_size) {
    int i,j,k,p,jfirst;
    Ullong key;
    Triel tt;
    
    for (i=0; i<npts; i++) {
        int ii;
        if (i>=vrtx_size) {
            ii = i % vrtx_size;
        }
        Vertex vrt(ii);
        vertices.push_back(vrt);
    }
    
    // Create a table so that given a point number, we can find one triangle with it as a vertex. Wrap a 3D triangle to a 1D point list.
    for (j=0; j<ntree; j++) {
        if (thelist[j].stat <= 0) continue;
        tt = thelist[j];
        for (k=0; k<3; k++) trindx[tt.p[k]] = j;
    }
    
    // Now loop over the sites
    for (p=0; p<npts; p++) {
        
        tt = thelist[trindx[p]];
        
        // Get the vertices into canonical order
        if (tt.p[0] == p) {i = tt.p[1]; j = tt.p[2];}
        else if (tt.p[1] == p) {i = tt.p[2]; j = tt.p[0];}
        else if (tt.p[2] == p) {i = tt.p[0]; j = tt.p[1];}
        else throw("triangle should contain p");
        
        // Save starting index
        jfirst = j;
        
        double dx = pvec[i].x[0] - pvec[p].x[0];
        double dy = pvec[i].x[1] - pvec[p].x[1];
        double dist = dx*dx + dy*dy;
        int pp = p % vrtx_size;
        if( dist < threshold ){
            int ii = i % vrtx_size;
            vertices[pp].addNeigh( ii );
        }
        
        dx = pvec[j].x[0] - pvec[p].x[0];
        dy = pvec[j].x[1] - pvec[p].x[1];
        dist = dx*dx + dy*dy;
        if( dist < threshold ){
            int jj = j % vrtx_size;
            vertices[pp].addNeigh( jj );
        }
        
        // Go around counterclockwise
        while (1) {
            key = hashfn.int64(i) - hashfn.int64(p);
            if ( ! linehash->get(key,k) ) throw("Delaunay is incomplete");
            
            double dx = pvec[k].x[0] - pvec[p].x[0];
            double dy = pvec[k].x[1] - pvec[p].x[1];
            double dist = dx*dx + dy*dy;
            if( dist < threshold ){
                int kk = k % vrtx_size;
                vertices[pp].addNeigh( kk );
            }
            
            if (k == jfirst) break;         // circumnavigation completed, normal way out
            j=i;
            i=k;
        }
    }
}


// Convex hull
struct Convexhull : Delaunay {
    int nhull;
    int *hullpts;
    Convexhull(vector<Point<2> > pvec);
};

Convexhull::Convexhull(vector<Point<2> > pvec) : Delaunay(pvec,2), nhull(0) {
    int i,j,k,pstart;
    vector<int> nextpt(npts);
    for (j=0; j<ntree; j++) {
        if (thelist[j].stat != -1) continue;
        for (i=0,k=1; i<3; i++,k++) {
            if (k == 3) k=0;
            if (thelist[j].p[i] < npts && thelist[j].p[k] < npts) break;
        }
        if (i==3) continue;
        ++nhull;
        nextpt[(pstart = thelist[j].p[k])] = thelist[j].p[i];
    }
    if (nhull == 0) throw("no hull segments found");
    hullpts = new int[nhull];
    j=0;
    i = hullpts[j++] = pstart;
    while ((i=nextpt[i]) != pstart) hullpts[j++] = i;
}



// T1 topological space
struct T1set {
    T1set() {}
    int n1, n2;
    int non1, non2;
    void addElem(int i, int j);
    void addPair(int i, int p);
};

void T1set::addElem(int i, int j) {
    n1 = i;
    n2 = j;
}

void T1set::addPair(int i, int p) {
    if (p == 1){
        non1 = i;
    }
    else {
        non2 = i;
    }
}




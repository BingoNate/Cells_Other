
// +++++++
// ++ Adapted from Numerical Recipes 3
// +++++++

#pragma once 

#include "pointbox.h"
#include "basic.h"
#include <vector>
#include <cmath>

using namespace std;


// Structure of the kd-tree
template<int DIM> struct KDtree {
    static const double BIG;                    // size of the root box, value set below
    int nboxes, npts;                           // number of boxes, number of points
    vector<Point<DIM> > &ptss;                  // keep a reference to the points making up the tree
    Boxnode<DIM> *boxes;                        // the array of Boxnodes that form the tree
    vector<int> ptindx, rptindx;                // index of points and reverse index
    double *coords;                             // point coordinates rearranged contiguously
    KDtree( vector<Point<DIM> > &pts );         // constructor
    ~KDtree() { delete [] boxes; }              // destructor
    
    // utility functions
    double disti( int jpt, int kpt );
    int locate( Point<DIM> pt );
    int locate( int jpt );
    
    // application functions
    int nearest( Point<DIM> pt );
    void nnearest( int jpt, int *nn, double *dn, int n );
    static void sift_down( double *heap, int *ndx, int nn );        // used by nnearest
    int locatenear( Point<DIM> pt, double r, int *list, int nmax );
};


template<int DIM> const double KDtree<DIM>::BIG(1.0e+99);


// Function to permute indx[0..n-1] to make arr[indx[0..k-1]] <= arr[indx[k]] <= arr[indx[k+1..n-1]]. The array arr is not modified. Used method is quickselect (partioning).
int selecti( const int k, int *indx, int n, double *arr ){
    
    int i_mid;
    double a;
    int i_first = 0;
    int i_last = n-1;
    
    for(;;){
        
        if( i_last <= i_first+1 ){                                                  // Active partition contains 1 or 2 elements
            if( i_last == i_first+1 && arr[indx[i_last]] < arr[indx[i_first]] )     // Active partition contains 2 elements
                SWAP( indx[i_first], indx[i_last] );
            return indx[k];
        }
        
        else{
            i_mid = (i_first+i_last)/2;
            SWAP( indx[i_mid], indx[i_first+1] );
            if( arr[indx[i_first]] > arr[indx[i_last]] ) SWAP( indx[i_first], indx[i_last] );
            if( arr[indx[i_first+1]] > arr[indx[i_last]] ) SWAP( indx[i_first+1], indx[i_last] );
            if( arr[indx[i_first]] > arr[indx[i_first+1]] ) SWAP( indx[i_first], indx[i_first+1] );
            int i = i_first+1;
            int j = i_last;
            int ia = indx[i_first+1];
            a = arr[ia];
            for(;;){
                do i++; while( arr[indx[i]] < a );
                do j--; while( arr[indx[j]] > a );
                if(j < i) break;
                SWAP( indx[i], indx[j] );
            }
            indx[i_first+1] = indx[j];
            indx[j] = ia;
            if( j >= k ) i_last = j-1;
            if( j <= k ) i_first = i;
            
        }
        
    }
    
}


// Construct a KD tree from a vector of points
template<int DIM> KDtree<DIM>::KDtree( vector<Point<DIM> > &pts ) : ptss(pts), npts(pts.size()), ptindx(npts), rptindx(npts) {
    
    int ntmp, m, k, kk, j, nowtask, jbox, np, tmom, tdim, ptlo, pthi;
    int *hp;
    double *cp;
    int taskmom[50], taskdim[50];
    for(k = 0; k < npts; k++) ptindx[k] = k;    // Initialize the index of points
    
    // Calculate the number of boxes and allocate memory for them
    m = 1;
    for( ntmp = npts; ntmp; ntmp /= 2 ){
        m *= 2;
    }
    nboxes = 2*npts - m/2;
    if( m < nboxes ) nboxes = m;
    nboxes--;
    boxes = new Boxnode<DIM>[nboxes];
    
    // Copy the coordinates into a contiguous array
    coords = new double[DIM*npts];
    for( j = 0, kk = 0; j < DIM; j++, kk += npts ){
        for( k = 0; k < npts; k++ ) coords[kk+k] = pts[k].x[j];
    }
    
    // Initialize the root box and put it on the task list for subdivision
    Point<DIM> lo(-BIG,-BIG,-BIG), hi(BIG,BIG,BIG);
    boxes[0] = Boxnode<DIM>(lo,hi,0,0,0,0,npts-1);
    jbox = 0;
    taskmom[1] = 0;       // which box
    taskdim[1] = 0;       // which dimension
    nowtask = 1;
    while (nowtask){
        tmom = taskmom[nowtask];
        tdim = taskdim[nowtask--];
        ptlo = boxes[tmom].ptlo;
        pthi = boxes[tmom].pthi;
        hp = &ptindx[ptlo];             // points to the left end of subdivision
        cp = &coords[tdim*npts];        // points to coordinate list for current dim.
        np = pthi - ptlo + 1;           // number of points in the subdivision
        kk = (np-1)/2;                  // index of last point on left (boundary point)
        (void) selecti(kk,hp,np,cp);    // here is where all the work is done!
        
        // Now create the daughters and push them onto the task list if they need further subdivision
        hi = boxes[tmom].hi;
        lo = boxes[tmom].lo;
        hi.x[tdim] = lo.x[tdim] = coords[tdim*npts + hp[kk]];
        boxes[++jbox] = Boxnode<DIM>(boxes[tmom].lo,hi,tmom,0,0,ptlo,ptlo+kk);
        boxes[++jbox] = Boxnode<DIM>(lo,boxes[tmom].hi,tmom,0,0,ptlo+kk+1,pthi);
        boxes[tmom].dau1 = jbox-1;
        boxes[tmom].dau2 = jbox;
        if (kk > 1){
            taskmom[++nowtask] = jbox-1;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
        if (np-kk > 3){
            taskmom[++nowtask] = jbox;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
    }
    
    for( j = 0; j < npts; j++ ) rptindx[ptindx[j]] = j;     // Create reverse index
    delete [] coords;
    
}


// Returns the distance between two points in the kdtree given their indices in the array of points, but returns a large value if the points are identical
template<int DIM> double KDtree<DIM>::disti(int jpt, int kpt){
    if (jpt == kpt) return BIG;
    else return dist(ptss[jpt], ptss[kpt]);
}


// Given an arbitrary point pt, return the index of which kdtree box it is in
template<int DIM> int KDtree<DIM>::locate(Point<DIM> pt){
    int nb,d1,jdim;
    nb = jdim = 0;                  // Start with the root box
    while (boxes[nb].dau1) {        // As far as possible down the tree
        d1 = boxes[nb].dau1;
        if (pt.x[jdim] <= boxes[d1].hi.x[jdim]) {nb = d1;}
        else {nb = boxes[nb].dau2;}
        jdim = ++jdim % DIM;        // Increment the dimension cyclically
    }
    return nb;
}


// Given the index of a point in the kdtree, return the index of which box it is in
template<int DIM> int KDtree<DIM>::locate(int jpt){
    int nb,d1,jh;
    jh = rptindx[jpt];      // The reverse index tells where the point lies in the index of points
    nb = 0;
    while (boxes[nb].dau1) {
        d1 = boxes[nb].dau1;
        if (jh <= boxes[d1].pthi) {nb = d1;}
        else {nb = boxes[nb].dau2;}
    }
    return nb;
}



// Given an arbitrary location pt, return the index of the nearest point in the kdtree
template<int DIM> int KDtree<DIM>::nearest(Point<DIM> pt){
    int i, k, nrst, ntask;
    int task[50];               // Stack for boxes waiting to be opened
    double dnrst = BIG, d;
    
    // First stage, we find the nearest kdtree point in same box as pt
    k = locate(pt);
    for(i=boxes[k].ptlo; i<=boxes[k].pthi; i++){        // Find nearest
        d = dist(ptss[ptindx[i]],pt);
        if(d < dnrst){
            nrst = ptindx[i];
            dnrst = d;
        }
    }
    
    // Second stage, we traverse the tree opening only possibly better boxes
    task[1] = 0;
    ntask = 1;
    while (ntask){
        k = task[ntask--];
        if( dist(boxes[k],pt) < dnrst ){    // Distance to closest point in box
            if( boxes[k].dau1 ){            // If not an end node, put on task list
                task[++ntask] = boxes[k].dau1;
                task[++ntask] = boxes[k].dau2;
            } else {
                for(i=boxes[k].ptlo; i<=boxes[k].pthi; i++){
                    d = dist(ptss[ptindx[i]],pt);
                    if(d < dnrst){
                        nrst = ptindx[i];
                        dnrst = d;
                    }
                }
            }
        }
    }
    return nrst;
}



// Given the index jpt of a point in a kdtree, return a list nn[0..n-1] of the indices of the n points in the tree nearest to point j, and a list dd[0..n-1] of their distances
template<int DIM> void KDtree<DIM>::nnearest(int jpt, int *nn, double *dn, int n){
    
    int i,k,ntask,kp;
    int task[50];       // Stack for boxes to be opened
    double d;
    if( n > npts-1 ) throw("Too many neighbors requested\n");
    for(i = 0; i < n; i++) dn[i] = BIG;
    
    // Find smallest mother box with enough points to initialize the heap
    kp = boxes[locate(jpt)].mom;
    while( boxes[kp].pthi - boxes[kp].ptlo < n ) kp = boxes[kp].mom;
    
    // Examine its points and save the n closest
    for( i = boxes[kp].ptlo; i <= boxes[kp].pthi; i++ ){
        
        if( jpt == ptindx[i] ) continue;
        d = disti(ptindx[i],jpt);
        if( d < dn[0] ){
            dn[0] = d;
            nn[0] = ptindx[i];
            if( n > 1 ) sift_down(dn,nn,n);     // Maintain the heap structure
        }
    }
    
    // Now traverse the tree opening only possibly better boxes
    task[1] = 0;
    ntask = 1;
    while( ntask ){
        k = task[ntask--];
        if( k == kp ) continue;     // Don't redo the box used to initialize
        if( dist(boxes[k],ptss[jpt]) < dn[0] ){
            if( boxes[k].dau1 ){
                task[++ntask] = boxes[k].dau1;
                task[++ntask] = boxes[k].dau2;
            } else{
                for( i = boxes[k].ptlo; i <= boxes[k].pthi; i++ ){
                    d = disti(ptindx[i],jpt);
                    if( d < dn[0] ){
                        dn[0] = d;
                        nn[0] = ptindx[i];
                        if( n > 1 ) sift_down(dn,nn,n);     // Maintain the heap
                    }
                }
            }
        }
        
    }
    
    return;
    
}



// Fix heap[0..nn-1] whose first element (only) may be wrongly filled. Make a corresponding permutation in ndx[0..nn-1].
template<int DIM> void KDtree<DIM>::sift_down(double *heap, int *ndx, int nn){
    
    int n = nn-1;
    int j,jold,ia;
    double a;
    a = heap[0];
    ia = ndx[0];
    jold = 0;
    j = 1;
    
    while( j <= n ){
        if( j < n && heap[j] < heap[j+1] ) j++;
        if( a >= heap[j] ) break;
        heap[jold] = heap[j];
        ndx[jold] = ndx[j];
        jold = j;
        j = 2*j + 1;
    }
    
    heap[jold] = a;
    ndx[jold] = ia;
    
}



// Given a point pt and radius r, returns a value nret such that list[0..nret-1] is a list of all kdtree points within a radius r of pt, up to a user-specified maximum of nmax points.
template<int DIM> int KDtree<DIM>::locatenear(Point<DIM> pt, double r, int *list, int nmax){
    
    int k,i,nb,nbold,nret,ntask,jdim,d1,d2;
    int task[50];
    nb = jdim = nret = 0;
    if( r < 0.0 ) throw("Radius must be nonnegative\n");
    
    // Find the smallest box that contains the "ball" of radius r
    while( boxes[nb].dau1 ){
        nbold = nb;
        d1 = boxes[nb].dau1;
        d2 = boxes[nb].dau2;
        
        // Only need to check the dimension that divides the daughters
        if( pt.x[jdim] + r <= boxes[d1].hi.x[jdim] ) nb = d1;
        else if( pt.x[jdim] - r >= boxes[d2].lo.x[jdim] ) nb = d2;
        jdim = ++jdim % DIM;
        if( nb == nbold ) break;        // Neither daughter encloses the ball
        
    }
    
    // Now traverse the tree below the starting box only as needed
    task[1] = nb;
    ntask = 1;
    while( ntask ){
        
        k = task[ntask--];
        if( dist(boxes[k],pt) > r ) continue;          // Box and ball are disjoint
        if( boxes[k].dau1 ){                           // Expand box further when possible
            task[++ntask] = boxes[k].dau1;
            task[++ntask] = boxes[k].dau2;
        } else {                                        // Otherwise process points in the box
            for( i = boxes[k].ptlo; i <= boxes[k].pthi; i++ ){
                if( dist(ptss[ptindx[i]],pt) <= r && nret < nmax )
                    list[nret++] = ptindx[i];
                if( nret == nmax ) return nmax;         // Not enough space!
            }
        }
        
    }
    
    return nret;
}










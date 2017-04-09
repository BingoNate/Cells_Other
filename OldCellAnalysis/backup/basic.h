
#pragma once
#include <vector>
#include <cmath>

inline int ifloor( double doubval ){
    return ( doubval >= 0. ? (int)(doubval) : ((int)doubval)-1 );
}

inline long double dnearbyint( double doubval ){
    return ( doubval >= 0. ? (long)(doubval + 0.5) : (long)(doubval - 0.5) );
}

inline int inearbyint( double doubval ){
    return ( doubval >= 0. ? (int)(doubval + 0.5) : (int)(doubval - 0.5) );
}

template<class T>
inline void SWAP(T &a, T &b){
    T dum=a;
    a=b;
    b=dum;
}
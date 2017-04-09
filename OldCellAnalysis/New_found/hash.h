
// +++++++
// ++ Adapted from Numerical Recipes 3
// +++++++

#pragma once

#include "pointbox.h"
#include "basic.h"
#include <vector>
#include <cmath>

using namespace std;


typedef unsigned long long int Ullong;    // 64 bit integer


// Random hash function
struct Ranhash {
    inline Ullong int64(Ullong u) {
        Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
        v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
        v *= 4768777513237032717LL;
        v ^= v << 20; v ^= v >> 41; v ^= v << 5;
        return  v;
    }
    inline unsigned int int32(Ullong u)
    { return (unsigned int)(int64(u) & 0xffffffff) ; }
    inline double doub(Ullong u)
    { return 5.42101086242752217E-20 * int64(u); }
};


// Structure of the hash table
template<class keyT, class hfnT> struct Hashtable {
    int nhash, nmax, nn, ng;
    vector<int> htable, next, garbg;
    vector<Ullong> thehash;
    hfnT hash;
    Hashtable(int nh, int nv);
    int iget(const keyT &key);
    int iset(const keyT &key);
    int ierase(const keyT &key);
    int ireserve();
    int irelinquish(int k);
};


// Constructor
template<class keyT, class hfnT>
Hashtable<keyT,hfnT>::Hashtable(int nh, int nv): hash(sizeof(keyT)), nhash(nh), nmax(nv), nn(0), ng(0), htable(nh), next(nv), garbg(nv), thehash(nv) {
    
        for (int j=0; j<nh; j++) { htable[j] = -1; }
    
}


template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::iget(const keyT &key) {
    
    int j,k;
    Ullong pp = hash.fn(&key);
    j = (int)(pp % nhash);
    
    for (k = htable[j]; k != -1; k = next[k]) {
        if (thehash[k] == pp) {
            return k;
        }
    }
    
    return -1;
}


template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::iset(const keyT &key) {
    
    int j,k,kprev;
    Ullong pp = hash.fn(&key);
    j = (int)(pp % nhash);
    
    if (htable[j] == -1) {
        k = ng ? garbg[--ng] : nn++ ;
        htable[j] = k;
    } else {
        for (k = htable[j]; k != -1; k = next[k]) {
            if (thehash[k] == pp) {
                return k;
            }
            kprev = k;
        }
        k = ng ? garbg[--ng] : nn++ ;
        next[kprev] = k;
    }
    
    if (k >= nmax) throw("storing too many values");
    
    thehash[k] = pp;
    next[k] = -1;
    return k;
}


template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::ierase(const keyT &key) {
    
    int j,k,kprev;
    Ullong pp = hash.fn(&key);
    j = (int)(pp % nhash);
    if (htable[j] == -1) return -1;
    kprev = -1;
    
    for (k = htable[j]; k != -1; k = next[k]) {
        if (thehash[k] == pp) {
            if (kprev == -1) htable[j] = next[k];
            else next[kprev] = next[k];
            garbg[ng++] = k;
            return k;
        }
        kprev = k;
    }
    
    return -1;
}


template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::ireserve() {
    int k = ng ? garbg[--ng] : nn++ ;
    if (k >= nmax) throw("reserving too many values");
    next[k] = -2;
    return k;
}


template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::irelinquish(int k) {
    if (next[k] != -2) {return -1;}
    garbg[ng++] = k;
    return k;
}


// Structure of the hash function
template<class keyT, class elT, class hfnT>
struct Hash : Hashtable<keyT, hfnT> {
    using Hashtable<keyT,hfnT>::iget;
    using Hashtable<keyT,hfnT>::iset;
    using Hashtable<keyT,hfnT>::ierase;
    vector<elT> els;
    
    Hash(int nh, int nm) : Hashtable<keyT, hfnT>(nh, nm), els(nm) {}
    
    void set(const keyT &key, const elT &el)
    {els[iset(key)] = el;}
    
    int get(const keyT &key, elT &el) {
        int ll = iget(key);
        if (ll < 0) return 0;
        el = els[ll];
        return 1;
    }
    
    elT& operator[] (const keyT &key) {
        int ll = iget(key);
        if (ll < 0) {
            ll = iset(key);
            els[ll] = elT();
        }
        return els[ll];
    }
    
    int count(const keyT &key) {
        int ll = iget(key);
        return (ll < 0 ? 0 : 1);
    }
    
    int erase(const keyT &key) {
        return (ierase(key) < 0 ? 0 : 1);
    }
};


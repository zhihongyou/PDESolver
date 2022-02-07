#ifndef DATATYPE_H
#define DATATYPE_H

#include <iostream> 
#include <vector>
#include <cstdio>
#include <string>

template <typename T>
struct vector2d {
    T x;
    T y;
};

template <typename T>
struct vector3d {
    T x;
    T y;
    T z;
};

template <typename T>
struct tensor2d {
    T xx;
    T xy;
    T yx;
    T yy;
};


#endif

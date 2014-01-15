#ifndef __TETVOL_H__
#define __TETVOL_H__

#include "mainhead.h"

//======================================
//Calculate volume of a tetraheron  
//======================================
__device__ float tetVol(float *r){

float x1 = r[0];
float y1 = r[1];
float z1 = r[2];
float x2 = r[3];
float y2 = r[4];
float z2 = r[5];
float x3 = r[6];
float y3 = r[7];
float z3 = r[8];
float x4 = r[9];
float y4 = r[10];
float z4 = r[11];

float a11=x1-x2;
float a12=y1-y2;
float a13=z1-z2;

float a21=x2-x3;
float a22=y2-y3;
float a23=z2-z3;

float a31=x3-x4;
float a32=y3-y4;
float a33=z3-z4;

float vol0=a11*a22*a33+a12*a23*a31+a13*a21*a32;
  vol0 = vol0-a13*a22*a31-a11*a23*a32-a12*a21*a33;
  vol0=vol0/6.0;

return abs(vol0);

}//tetVol



#endif 

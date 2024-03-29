#ifndef __FORCECALC_H__
#define __FORCECALC_H__

#include "mainhead.h"
#include "parameters.h"
#include "UserDefined.h"
#include "invert4x4.h"
#include "setMetric.h"


//=============================================================
//calculate forces on all 4 nodes in tetrahedra
//=============================================================
__device__ void force_calc(float *r0,float *r,float (&F)[12],int *TetNodeRank,float *pe,int mytet,float myVol,float t,float *swell,int tetID){

	float u[4],v[4],w[4];
	float eps[9];
	float a[4] = {0.0};
	float b[4] = {0.0};
	float c[4] = {0.0};
	float localPe = 0.0;
	float A[16] = {0.0};
  float Ainv[16] = {0.0};
  float stress[9];
  float eps_trace;

  //do swelling
  sweller(r0,swell,tetID);

	//clacluate displacements from original position and zero out forces
	for(int n=0;n<4;n++){
    A[n*4] = 1.0;
    A[n*4+1] = r0[n*3];
    A[n*4+2] = r0[n*3+1];
    A[n*4+3] = r0[n*3+2];
		u[n] = (r[n*3]-r0[n*3]);
		v[n] = (r[1+n*3]-r0[1+n*3]);
		w[n] = (r[2+n*3]-r0[2+n*3]);
		F[0+n*3] = 0.0;
		F[1+n*3] = 0.0;
		F[2+n*3] = 0.0;
	}//n

  // invert A
  gluInvertMatrix(A,Ainv);


	//matrix multipy Ainv and u,v,w to get shape funcitons
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			a[i]+=Ainv[4*i+j]*u[j];
			b[i]+=Ainv[4*i+j]*v[j];
			c[i]+=Ainv[4*i+j]*w[j];
		}//j
	}//i

	//noe we calculate epsilon tensor
	eps[3*0+0] = a[1]+0.5*(a[1]*a[1]+b[1]*b[1]+c[1]*c[1]);
	eps[3*1+1] = b[2]+0.5*(a[2]*a[2]+b[2]*b[2]+c[2]*c[2]);
	eps[3*2+2] = c[3]+0.5*(a[3]*a[3]+b[3]*b[3]+c[3]*c[3]);
	eps[3*0+1] = 0.5*(a[2]+b[1]+a[1]*a[2]+b[1]*b[2]+c[1]*c[2]);
	eps[3*1+0] = eps[3*0+1];
	eps[3*0+2] = 0.5*(a[3]+c[1]+a[1]*a[3]+b[1]*b[3]+c[1]*c[3]);
	eps[3*2+0] = eps[3*0+2];
	eps[3*1+2] = 0.5*(b[3]+c[2]+a[2]*a[3]+b[2]*b[3]+c[2]*c[3]);
	eps[3*2+1] = eps[3*1+2];

  //calculate trace of epsilon
  eps_trace  = eps[3*0+0]+eps[3*1+1]+eps[3*2+2];

  //calculate stress tensor
  stress[3*0+0] = LAMBDA*eps_trace+2.0*MU*eps[3*0+0];
  stress[3*0+1] = 2.0*MU*eps[3*0+1];
  stress[3*0+2] = 2.0*MU*eps[3*0+2];
  stress[3*1+0] = stress[3*0+1];
  stress[3*1+1] = LAMBDA*eps_trace+2.0*MU*eps[3*1+1];
  stress[3*1+2] = 2.0*MU*eps[3*1+2];
  stress[3*2+0] = stress[3*0+2]; 
  stress[3*2+1] = stress[3*1+2];
  stress[3*2+2] = LAMBDA*eps_trace+2.0*MU*eps[3*2+2];
	
  //update swelling 
  if(t>0.1){ updateSwell(swell,stress,tetID);}

  //calculate potential energy
	localPe += cxxxx*(eps[3*0+0]*eps[3*0+0]+eps[3*1+1]*eps[3*1+1]+eps[3*2+2]*eps[3*2+2]);
	localPe += 2.0*cxxyy*(eps[3*0+0]*eps[3*1+1]+eps[3*1+1]*eps[3*2+2]+eps[3*0+0]*eps[3*2+2]);
	localPe += 4.0*cxyxy*(eps[3*0+1]*eps[3*0+1]+eps[3*1+2]*eps[3*1+2]+eps[3*2+0]*eps[3*2+0]);

		//send potential to global memory
	pe[mytet] = localPe*myVol;

	//now can calculate forces
	for(int n = 0;n<4;n++){
		//x
		
		F[0+n*3]+=-cxxxx*2.0*eps[3*0+0]*Ainv[4*1+n]*(1.0+a[1]);
		F[0+n*3]+=-cxxxx*2.0*eps[3*1+1]*Ainv[4*2+n]*a[2];
		F[0+n*3]+=-cxxxx*2.0*eps[3*2+2]*Ainv[4*3+n]*a[3];
		F[0+n*3]+=-2.0*cxxyy*eps[3*0+0]*Ainv[4*2+n]*a[2];
		F[0+n*3]+=-2.0*cxxyy*eps[3*1+1]*Ainv[4*1+n]*(1.0+a[1]);
		F[0+n*3]+=-2.0*cxxyy*eps[3*1+1]*Ainv[4*3+n]*a[3];
		F[0+n*3]+=-2.0*cxxyy*eps[3*2+2]*Ainv[4*2+n]*a[2];
		F[0+n*3]+=-2.0*cxxyy*eps[3*2+2]*Ainv[4*1+n]*(1.0+a[1]);
		F[0+n*3]+=-2.0*cxxyy*eps[3*0+0]*Ainv[4*3+n]*a[3];
		F[0+n*3]+=-4.0*cxyxy*eps[3*0+1]*Ainv[4*2+n]*(1.0+a[1]);
		F[0+n*3]+=-4.0*cxyxy*eps[3*0+1]*Ainv[4*1+n]*a[2];
		F[0+n*3]+=-4.0*cxyxy*eps[3*1+2]*Ainv[4*3+n]*a[2];
		F[0+n*3]+=-4.0*cxyxy*eps[3*1+2]*Ainv[4*2+n]*a[3];
		F[0+n*3]+=-4.0*cxyxy*eps[3*2+0]*Ainv[4*3+n]*(1.0+a[1]);
		F[0+n*3]+=-4.0*cxyxy*eps[3*2+0]*Ainv[4*1+n]*a[3];

		//y
		
		F[1+n*3]+=-cxxxx*2.0*eps[3*0+0]*Ainv[4*1+n]*b[1];
		F[1+n*3]+=-cxxxx*2.0*eps[3*1+1]*Ainv[4*2+n]*(1.0+b[2]);
		F[1+n*3]+=-cxxxx*2.0*eps[3*2+2]*Ainv[4*3+n]*b[3];
		F[1+n*3]+=-2.0*cxxyy*eps[3*0+0]*Ainv[4*2+n]*(1.0+b[2]);
		F[1+n*3]+=-2.0*cxxyy*eps[3*1+1]*Ainv[4*1+n]*b[1];
		F[1+n*3]+=-2.0*cxxyy*eps[3*1+1]*Ainv[4*3+n]*b[3];
		F[1+n*3]+=-2.0*cxxyy*eps[3*2+2]*Ainv[4*2+n]*(1.0+b[2]);
		F[1+n*3]+=-2.0*cxxyy*eps[3*2+2]*Ainv[4*1+n]*b[1];
		F[1+n*3]+=-2.0*cxxyy*eps[3*0+0]*Ainv[4*3+n]*b[3];
		F[1+n*3]+=-4.0*cxyxy*eps[3*0+1]*Ainv[4*1+n]*(1.0+b[2]);
		F[1+n*3]+=-4.0*cxyxy*eps[3*0+1]*Ainv[4*2+n]*b[1];
		F[1+n*3]+=-4.0*cxyxy*eps[3*1+2]*Ainv[4*3+n]*(1.0+b[2]);
		F[1+n*3]+=-4.0*cxyxy*eps[3*1+2]*Ainv[4*2+n]*b[3];
		F[1+n*3]+=-4.0*cxyxy*eps[3*2+0]*Ainv[4*3+n]*b[1];
		F[1+n*3]+=-4.0*cxyxy*eps[3*2+0]*Ainv[4*1+n]*b[3];

		//z
		
		F[2+n*3]+=-cxxxx*2.0*eps[3*0+0]*Ainv[4*1+n]*c[1];
		F[2+n*3]+=-cxxxx*2.0*eps[3*1+1]*Ainv[4*2+n]*c[2];
		F[2+n*3]+=-cxxxx*2.0*eps[3*2+2]*Ainv[4*3+n]*(1.0+c[3]);
		F[2+n*3]+=-2.0*cxxyy*eps[3*0+0]*Ainv[4*2+n]*c[2];
		F[2+n*3]+=-2.0*cxxyy*eps[3*1+1]*Ainv[4*1+n]*c[1];
		F[2+n*3]+=-2.0*cxxyy*eps[3*1+1]*Ainv[4*3+n]*(1.0+c[3]);
		F[2+n*3]+=-2.0*cxxyy*eps[3*2+2]*Ainv[4*2+n]*c[2];
		F[2+n*3]+=-2.0*cxxyy*eps[3*2+2]*Ainv[4*1+n]*c[1];
		F[2+n*3]+=-2.0*cxxyy*eps[3*0+0]*Ainv[4*3+n]*(1.0+c[3]);
		F[2+n*3]+=-4.0*cxyxy*eps[3*0+1]*Ainv[4*2+n]*c[1];
		F[2+n*3]+=-4.0*cxyxy*eps[3*0+1]*Ainv[4*1+n]*c[2];
		F[2+n*3]+=-4.0*cxyxy*eps[3*1+2]*Ainv[4*2+n]*(1.0+c[3]);
		F[2+n*3]+=-4.0*cxyxy*eps[3*1+2]*Ainv[4*3+n]*c[2];
		F[2+n*3]+=-4.0*cxyxy*eps[3*2+0]*Ainv[4*1+n]*(1.0+c[3]);
		F[2+n*3]+=-4.0*cxyxy*eps[3*2+0]*Ainv[4*3+n]*c[1];

		
		}//n




}//force_calc


//=============================================================
//calculate confining forces
//=============================================================
__device__ void confineForce(float *r,float (&F)[12]){
  int nLoc;
  for (int n=0;n<4;n++){
    nLoc = 2+n*3;
    if(r[nLoc]>CONFTOP){
      F[nLoc]+=-CONFFORCE; 
    }
    if(r[nLoc]<CONFBOTTOM){
      F[nLoc]+=CONFFORCE; 
    }//if
  }//for n
}


#endif//__FORCECALC_H__

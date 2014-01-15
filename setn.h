#ifndef __SETN_H__
#define __SETN_H__

#include "mainhead.h"
#include "UserDefined.h"
#include "initDirector.h"
#include "printDirector.h"
//#include "directorRelax.h"
#include <math.h>




//program to read in top and bottom boundary conditions
//from the top and bottom data files. CRUDE for nowi
//
/*
void getNfromFile(float (&botBC)[inX][inY][2],float (&topBC)[inX][inY][2]){


	int topI,topJ,botI,botJ;
	float topTheta,topPhi,botTheta,botPhi;
	FILE* botIN;
	FILE* topIN;
	
    botIN = fopen("Bounds//bottom.dat","r+");
	topIN = fopen("Bounds//top.dat","r+");
	
	for(int i=0;i<inX;i++){
		for(int j=0;j<inY;j++){
			fscanf(topIN,"%d %d %f %f\n",&topI,&topJ,&topTheta,&topPhi);
			fscanf(botIN,"%d %d %f %f\n",&botI,&botJ,&botTheta,&botPhi);
			topBC[topI-1][topJ-1][0] = topTheta;
			topBC[topI-1][topJ-1][1] = topPhi;
			botBC[botI-1][botJ-1][0] = botTheta;
			botBC[botI-1][botJ-1][1] = botPhi;
		}//j
	}//i
	

}//getNfromFile
*/

//set the theta and phi corresponding to the average director inside each tetrahedra
//simple now but will connect this to the mathematica GUI
//****must execute get_tet_pos first*****
void set_n(TetArray &i_Tet,int Ntets){
	float rx,ry,rz,theta=0.0,phi=0.0;
/*
	//zero out top and botom directo profiles
	float botBC[inX][inY][2] = {0.0};
	float topBC[inX][inY][2] = {0.0};


	//max and min values of mesh to determine scalling for 
	//director simulation
	float maxX = i_Tet.max(0), maxY = i_Tet.max(1), maxZ = i_Tet.max(2);
	float minX = i_Tet.min(0), minY = i_Tet.min(1), minZ = i_Tet.min(2);

	//get aspect ratio of mesh
	float xDim = maxX-minX, yDim = maxY-minY, zDim = maxZ-minZ;

	//get optimal Z dimension to have square lattice inX x inY x inZ
	//of the same aspect ratio of mesh
	int inZ = int(float(inX)*zDim/xDim);

	//create director array of correct size to simulate 
	float * THETA;
	float * PHI;
	THETA = (float *) malloc(inX*inY*inZ*sizeof(float));
	PHI = (float *) malloc(inX*inY*inZ*sizeof(float));

	//read in top and bottom director profiles
	getNfromFile(botBC,topBC);

	//initialize THETA and PHI to be sent to simulation using the 
	//boundary conditions read in from the file
	initN(THETA,PHI,inZ,topBC,botBC);

	//print director before sim
	printDirector(THETA,PHI,inZ,1);

	//simulate director relaxation
	directorRelax(THETA,PHI,inX,inY,inZ);

	//print director after sim
	printDirector(THETA,PHI,inZ,2);
*/
	
	for(int i=0;i<Ntets;i++){

		//get position of tetrahedra
		rx = i_Tet.get_pos(i,0);
		ry = i_Tet.get_pos(i,1);
		rz = i_Tet.get_pos(i,2);
		
		//turn positions into director
	//	getThPhi(rx,ry,rz,theta,phi,THETA,PHI,inZ,maxX,minX,maxY,minY,maxZ,minZ);

		//assign theta and phi to tetrahedra
		i_Tet.set_theta(i,theta);
		i_Tet.set_phi(i,phi);
	}//i


		

}//set angle


#endif//__SETN_H__

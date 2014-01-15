#ifndef __CONNECTIONGRAPH_H__
#define __CONNECTIONGRAPH_H__

#include <math.h>
#include "mainhead.h"
#include "setMetric.h"

//holds connectvity graph of surface nodes
//calculates gaussian curvature from data
class connectionGraph{
  public:
    int *nabs;
    float *position;
    float *initialPosition;
    int totalNodes;
    int maxNabs;
    int *surf;
    int *triangle;
    int totalTriangles;
    float *triangleArea;
    float *myTriangle;
    int maxTrianglesPerNode;
    int *nMyTriangle;
    int *reduceMap;
    float *gauss;
    float *turningAngle;
    float *metric;
    connectionGraph(int totalNodes);
   ~connectionGraph();
  
  void setPosition(int i, const float &x, const float &y, const float &z);
  void makeSurf(int i);
  void makeNeighbors(int n1, int n2);
  void print(const char * filename);
  void printFrame(int step);
  void addTriangle(const int n1, const int n2, const int n3);
  int amISurf(int n);
  void reduce();
  void calcGauss();
  void update(DevDataBlock *dev_dat, HostDataBlock *host_dat);
  void setInitial();
  void dumpState(int step);
};



//initialize suface connectivity graph
connectionGraph::connectionGraph(int N){
  totalNodes = N;
  maxNabs = 50;
  maxTrianglesPerNode = 40;
  totalTriangles = 0;
  position = new float[totalNodes*3];
  initialPosition = new float[totalNodes*3];
  nabs = new int[totalNodes*maxNabs];
  surf = new int[totalNodes];
  triangle = new int[totalNodes*3];
  triangleArea = new float[totalNodes];
  myTriangle = new float[totalNodes*maxTrianglesPerNode];
  nMyTriangle = new int[totalNodes];
  reduceMap = new int[totalNodes];
  gauss = new float[totalNodes];
  turningAngle = new float[totalNodes];
  metric = new float[totalNodes];

  for(int i=0;i<totalNodes*maxNabs;i++){
    nabs[i] = -10;
  }//i
  for(int j=0;j<totalNodes;j++){
    surf[j] = 0;
    triangleArea[j] = 0.0;
    nMyTriangle[j] = 0;
    reduceMap[j] = -10;
    gauss[j] = 0.0;
    turningAngle[j] = 0.0;
  }//j
  for(int k=0;k<totalNodes*3;k++){
    triangle[k] = -10;
  }//k
  for(int l=0;l<totalNodes*maxTrianglesPerNode;l++){
    myTriangle[l] = -10;
  }//l
}//connectionGraph::connectionGraph(int N){


//deconstructor
connectionGraph::~connectionGraph(){
  delete nabs;
  nabs = NULL;
  delete position;
  position = NULL;
  delete surf;
  surf = NULL;
}//~

//set initial state
void connectionGraph::setInitial(){
  float thetaSum,areaSum,theta,cosTheta;
  float Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,LB,LC,nx,ny,nz;
  int n1,n2,n3,t; 
  for(int n=0;n<totalNodes;n++){
    if(surf[n]){
      //set initial postions
      nx = position[n*3];
      ny = position[n*3+1];
      nz = position[n*3+2];
      initialPosition[n*3] = nx;
      initialPosition[n*3+1] = ny;
      initialPosition[n*3+2] = nz;
      metric[n] = getMetricLocal(nx,ny,nz);
        
      //calculate turnign angle
      thetaSum = 0.0;
      areaSum = 0.0;
      for(int nt=0;nt<nMyTriangle[n];nt++){
          t = myTriangle[n*maxTrianglesPerNode+nt];
          areaSum += triangleArea[t];
          n1 = triangle[t*3];    
          n2 = triangle[t*3+1];    
          n3 = triangle[t*3+2];
          if(n2==n){
            n2 = n1;
            n1 = n;
          }else if(n3 == n){
            n3 = n1;
            n1 = n;
          }//if || elsee
          Ax = position[n1*3];
          Ay = position[n1*3+1];
          Az = position[n1*3+2];
          Bx = position[n2*3]-Ax;
          By = position[n2*3+1]-Ay;
          Bz = position[n2*3+2]-Az;
          Cx = position[n3*3]-Ax;
          Cy = position[n3*3+1]-Ay;
          Cz = position[n3*3+2]-Az;
          
          LB = sqrt(Bx*Bx + By*By + Bz*Bz);
          LC = sqrt(Cx*Cx + Cy*Cy + Cz*Cz);
          
          cosTheta = (Bx*Cx+By*Cy+Bz*Cz)/(LB*LC);
          theta = acos(cosTheta);
          thetaSum+=theta; 
      }//n
      turningAngle[n] = thetaSum;
    }//if surf
  }//i

}//setTurning



//update positions from GPU
void connectionGraph::update(DevDataBlock *dev_dat
                            , HostDataBlock *host_dat){

  //for pitched memory transfor
  size_t height3 = 3;
  size_t widthNODE = totalNodes;

  HANDLE_ERROR( cudaMemcpy2D( host_dat->host_r
                          , widthNODE*sizeof(float)
                          , dev_dat->dev_r
                          , dev_dat->dev_rpitch
                          , widthNODE*sizeof(float)
                          , height3
                          , cudaMemcpyDeviceToHost ) );

  for(int n=0;n<totalNodes;n++){
    position[n*3] = host_dat->host_r[n+0*totalNodes];
    position[n*3+1] = host_dat->host_r[n+1*totalNodes];
    position[n*3+2] = host_dat->host_r[n+2*totalNodes];
  }//i
  
}//update

//calculate Gaussian curvature of every point
void connectionGraph::calcGauss(){
  float thetaSum,areaSum,theta,cosTheta;
  float Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,LB,LC;
  int n1,n2,n3,t; 
  float pi = 3.14159265359;
  for(int n=0;n<totalNodes;n++){
    if(surf[n]){
      thetaSum = 0.0;
      areaSum = 0.0;
      for(int nt=0;nt<nMyTriangle[n];nt++){
          t = myTriangle[n*maxTrianglesPerNode+nt];
          areaSum += triangleArea[t];
          n1 = triangle[t*3];    
          n2 = triangle[t*3+1];    
          n3 = triangle[t*3+2];
          if(n2==n){
            n2 = n1;
            n1 = n;
          }else if(n3 == n){
            n3 = n1;
            n1 = n;
          }//if || elsee
          Ax = position[n1*3];
          Ay = position[n1*3+1];
          Az = position[n1*3+2];
          Bx = position[n2*3]-Ax;
          By = position[n2*3+1]-Ay;
          Bz = position[n2*3+2]-Az;
          Cx = position[n3*3]-Ax;
          Cy = position[n3*3+1]-Ay;
          Cz = position[n3*3+2]-Az;
          
          LB = sqrt(Bx*Bx + By*By + Bz*Bz);
          LC = sqrt(Cx*Cx + Cy*Cy + Cz*Cz);
          
          cosTheta = (Bx*Cx+By*Cy+Bz*Cz)/(LB*LC);
          theta = acos(cosTheta);
          thetaSum+=theta; 
      }//n
      gauss[n] = 3.0*(turningAngle[n]-thetaSum)/areaSum;
    }//if surf
  }//i
}//calcgauss


void connectionGraph::reduce(){
  int count = 0;
  for(int i=0;i<totalNodes;i++){
    if(surf[i]){
      reduceMap[i] = count;
      count++;
    }//if surf
  }//i
}//reduce


void connectionGraph::addTriangle(const int n1, const int n2, const int n3){
  triangle[totalTriangles*3] = n1;
  triangle[totalTriangles*3+1] = n2;
  triangle[totalTriangles*3+2] = n3;

  //calculate area 
  float area = 0.0;
  float Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
  int mt1,mt2,mt3;
  Ax = position[n1*3];
  Ay = position[n1*3+1];
  Az = position[n1*3+2];
  Bx = position[n2*3];
  By = position[n2*3+1];
  Bz = position[n2*3+2];
  Cx = position[n3*3];
  Cy = position[n3*3+1];
  Cz = position[n3*3+2];
  
  area = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))*0.5;
  triangleArea[totalTriangles] = area;
  
  //add to triangle list for each node
  mt1 = n1*maxTrianglesPerNode+nMyTriangle[n1];
  mt2 = n2*maxTrianglesPerNode+nMyTriangle[n2];
  mt3 = n3*maxTrianglesPerNode+nMyTriangle[n3];
  myTriangle[mt1] = totalTriangles;
  nMyTriangle[n1]++;
  myTriangle[mt2] = totalTriangles;
  nMyTriangle[n2]++;
  myTriangle[mt3] = totalTriangles;
  nMyTriangle[n3]++;
 
  
  //increment total triangles
  totalTriangles++;
}//addTriangle

void connectionGraph::setPosition(int i, const float &x, const float &y, const float &z){
  position[i*3] = x;
  position[i*3+1] = y;
  position[i*3+2] = z;
}//connectionGraph::setPosition

int connectionGraph::amISurf(int n){
  return surf[n];
}//amISurf

void connectionGraph::makeSurf(int i){
   surf[i]=1; 
}//connectionGraph::makeSurf(


void connectionGraph::makeNeighbors(int n1, int n2){
  int nab1,nab2;
  if((surf[n1]+surf[n2])==2){
    for(int i = 0;i<maxNabs;i++){
      nab1 = nabs[n1*maxNabs+i];
      if(nab1==n2){
        i = maxNabs;
      }else if(nab1==-10){
        nabs[n1*maxNabs+i] = n2;
        i = maxNabs;
      }//
    } //for i
    for(int j = 0;j<maxNabs;j++){
      nab2 = nabs[n2*maxNabs+j];
      if(nab2==n1){
        j = maxNabs;
      }else if(nab2==-10){
        nabs[n2*maxNabs+j] = n1;
        j = maxNabs;
      }//
     } //for j
  }//if surf
}//connectionGraph::makeNeighbors

void connectionGraph::printFrame(int step){
  char fout[60];
  sprintf(fout,"VTKOUT2//mesh%d.vtk",step);
  print(fout);
}//printFrame


//dump all data to be analyzed
void connectionGraph::dumpState(int step){
  FILE*out;
  char fout[60];
  sprintf(fout,"STATEDATA//state%d.dat",step);
  out = fopen(fout,"w");

  for(int n=0;n<totalNodes;n++){
    if(surf[n]){
      fprintf(out,"%f %f %f %f %f %f %f %f\n"
          ,initialPosition[n*3]
          ,initialPosition[n*3+1]
          ,initialPosition[n*3+2]
          ,position[n*3]
          ,position[n*3+1]
          ,position[n*3+2]
          ,metric[n]
          ,gauss[n]); 
    }//if surf  
  }//n 
  fclose(out); 
}//dump state


void connectionGraph::print(const char *filename){
  //calculate surface points
  int surfacePoints = 0;
  for(int i=0;i<totalNodes;i++){
    surfacePoints += surf[i];
  }//for i
 
  FILE *out; 
  out = fopen(filename,"w");
  
  //header info
  fprintf(out,"# vtk DataFile Version 3.1\n");
  fprintf(out,"Surface Graph\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(out,"\n");
  fprintf(out,"POINTS %d FLOAT\n",surfacePoints);
  
  //positions
  for(int n=0;n<totalNodes;n++){
    if(surf[n]){
      fprintf(out,"%f %f %f\n",position[n*3],position[n*3+1],position[n*3+2]);
    }//if surf
  }//n
  fprintf(out,"\n");
  
  //cells
  fprintf(out,"CELLS %d %d\n",totalTriangles,4*totalTriangles);
  fprintf(out,"\n");
  for(int n=0;n<totalTriangles;n++){
    fprintf(out,"%d %d %d %d\n",3,reduceMap[triangle[n*3]],reduceMap[triangle[n*3+1]],reduceMap[triangle[n*3+2]]);
  }//for n
  fprintf(out,"\n");
  
  //cell types
  fprintf(out,"CELL_TYPES %d \n",totalTriangles);
  for(int n=0;n<totalTriangles;n++){
    fprintf(out,"5\n");
  }//n
  fprintf(out,"\n");

  //gaussian curvature at each point
  fprintf(out,"POINT_DATA %d\n",surfacePoints);
  fprintf(out,"SCALARS gaussianCurvature FLOAT 1\n");
  fprintf(out,"LOOKUP_TABLE default\n");
  for(int n=0;n<totalNodes;n++){
    if(surf[n]){
      fprintf(out,"%f\n",gauss[n]);
    }//if surf
  }//n
  fprintf(out,"\n");
  //close file
  fclose(out);
}//void connectionGraph::print

#endif// __CONNECTIONGRAPH_H__

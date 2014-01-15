#ifndef __SETMETRIC_H__
#define __SETMETRIC_H__

//#include "mainhead.h"

//===================================//
//set metric as function of position //
//===================================//
__device__ float getMetric(float rx,float ry,float rz){
  
 //flat sheet --sphere
  float cx = 0.0, cy = 0.0;
  float bigR = 15.0;
  float cSwell = 1.4;
  float dr2,dx,dy,temp,bigR2,metric;
  
  bigR2 = bigR*bigR;
  dx = cx-rx;
  dy = cy-ry;
  dr2 = dx*dx+dy*dy;
  temp = 1.0+dr2/bigR2; 
  metric = cSwell/temp;

return metric;
}//getMetric



//===================================//
// swell "equilibrium" state         //
//===================================//
__device__ void swell(float *r0,float t){
  float xmid=0.0,ymid=0.0,zmid=0.0;
  float metricQR; //cube root of metric
  float rSwell[12];  
  float swellTime = 0.030;
  float slowSwell;
  //find initial center of tetrahedron to be used for metric(r0)
  for(int n=0;n<4;n++){
    xmid += r0[3*n];
    ymid += r0[3*n+1];
    zmid += r0[3*n+2];
  }//n

  xmid = xmid/4.0;
  ymid = ymid/4.0;
  zmid = zmid/4.0;

  //calculate tetrahedron with one vertix @ origin
  for(int n=0;n<4;n++){
    rSwell[3*n] = r0[3*n]-r0[0]; 
    rSwell[3*n+1] = r0[3*n+1]-r0[1]; 
    rSwell[3*n+2] = r0[3*n+2]-r0[2]; 
  }//n

  //get swelling factor
  metricQR = getMetric(xmid,ymid,zmid);
  //metricQR = 1.3;
  
  //slowly swell 
  slowSwell = 1.0;
  if(t<swellTime){slowSwell = t/swellTime;}
  metricQR = 1.0+(metricQR-1.0)*slowSwell;

  //swell 
  for(int n=0;n<4;n++){
    r0[3*n] = r0[0]+metricQR*rSwell[3*n];
    r0[3*n+1] = r0[1]+metricQR*rSwell[3*n+1];
    r0[3*n+2] = r0[2]+metricQR*rSwell[3*n+2];
  }//n

}//swell


//===================================//
//set metric as function of position //
//===================================//
float getMetricLocal(float rx,float ry,float rz){
   
  
 //flat sheet --sphere
  float cx = 0.0, cy = 0.0;
  float bigR = 15.0;
  float cSwell = 1.4;
  float dr2,dx,dy,temp,bigR2,metric;
  
  bigR2 = bigR*bigR;
  dx = cx-rx;
  dy = cy-ry;
  dr2 = dx*dx+dy*dy;
  temp = 1.0+dr2/bigR2;  
  metric = cSwell/temp;

return metric;
}//getMetric


//========================================//
//output metric data                      //
//========================================//
void metricData(){
  float xMin = -15.0,xMax = 15.0,dx = 1.0;
  float yMin = -10.0,yMax = 10.0,dy = 1.0;
  FILE *out;
  out = fopen("Output//metric.dat","w");
  float xcur = xMin, ycur = yMin;
  while(xcur<xMax){
    ycur = yMin;
    while(ycur<yMax){
      fprintf(out,"%f %f %f\n",xcur,ycur,getMetricLocal(xcur,ycur,0.0));
      ycur +=dy;
    }//ycur<yMax
    xcur += dx;
  }//xcur<xmax
  fclose(out);
}//metric data

#endif 

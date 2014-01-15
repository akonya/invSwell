#ifndef __BUILDCONNECTIONGRAPH_H__
#define __BUILDCONNECTIONGRAPH_H__

#include "parameters.h"
#include <math.h>

//build connection graph and suface
void buildConnectionGraph(NodeArray &i_Node,TetArray &i_Tet, connectionGraph &graph,int Nnodes, int Ntets){
  //find position of top/bottom surface
  float zTop = -100.0;
  float zBottom = 100.0;
  float locz,dzTop,dzBottom,dz;
  int nbs[4];
  int triCount,nb;
  for(int i=0;i<Nnodes;i++){
    locz = i_Node.get_pos(i,2);
    if(locz>zTop){zTop = locz;}
    if(locz<zBottom){zBottom = locz;}
  }//i
  
  for(int i=0;i<Nnodes;i++){
    locz = i_Node.get_pos(i,2);
    dzTop = zTop-locz;
    dzBottom = zBottom-locz;
    //only do top surface for now
    if(abs(dzTop)<0.001){
      graph.makeSurf(i);
      graph.setPosition(i,i_Node.get_pos(i,0),i_Node.get_pos(i,1),i_Node.get_pos(i,2));
    }//if abs(dzTop)<
  }//i
  
  for(int i=0;i<Ntets;i++){
    triCount = 0;
    for(int j=0;j<4;j++){
       nb = i_Tet.get_nab(i,j);
       if(graph.amISurf(nb)){
          nbs[triCount]=nb;
          triCount++;
       }//if surf 
    }//j
    
    if(triCount==3){
      graph.addTriangle(nbs[0],nbs[1],nbs[2]);
    }//if triCount==3
/*
    graph.makeNeighbors(nbs[0],nbs[1]);
    graph.makeNeighbors(nbs[0],nbs[2]);
    graph.makeNeighbors(nbs[0],nbs[3]);
    graph.makeNeighbors(nbs[1],nbs[2]);
    graph.makeNeighbors(nbs[1],nbs[3]);
    graph.makeNeighbors(nbs[2],nbs[3]);
 */   
  }//i
}//  buildConnectionGraph

#endif // __BUILDCONNECTIONGRAPH_H__

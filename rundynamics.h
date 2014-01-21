#ifndef __RUNDYNAMICS_H__
#define __RUNDYNAMICS_H__
#include "mainhead.h"
#include "parameters.h"
#include "printVTKframe.h"
#include "gpuForce.h"
#include "updateKernel.h"
#include "getEnergy.h"
#include "setMetric.h"
#include "connectionGraph.h"


//This funciton handles all dynamics which will be run
void run_dynamics(DevDataBlock *data
					, HostDataBlock *host_data
          , connectionGraph *graph
					, int Ntets
					, int Nnodes
					,int *Syncin
					,int *Syncout
					,int *g_mutex){

	//==============================================================
	//file to write energies to
	//===============================================================
		FILE*Eout;
		Eout = fopen("Output\\EvsT.dat","w");
		float pE,kE;

  //===============================================================
  //output metric to a file 
  //===============================================================
  metricData(); 


	
	//=================================================================
	//claclulate number of blocks to be executed
	//=================================================================
	cudaDeviceProp dev_prop;
	HANDLE_ERROR(cudaGetDeviceProperties(&dev_prop,0));
	int Threads_Per_Block = TPB;
	int BlocksTet = (Ntets+Threads_Per_Block)/Threads_Per_Block;
	int BlocksNode = (Nnodes+Threads_Per_Block)/Threads_Per_Block;

	printf("execute dynamnics kernel using:\n%d blocks\n%d threads per bock\n",BlocksTet,Threads_Per_Block);
	



	size_t widthTETS = Ntets;
	size_t height16 = 16;

	//================================================================
	// create start and stop events to measure performance
	//================================================================
	
	cudaEvent_t start, stop; 
	float elapsedTime;
	

	//================================================================
	// Begin Dynamics
	//================================================================

	for(int iKern=0;iKern<NSTEPS;iKern++){
	//time each kernal launch
	HANDLE_ERROR(cudaEventCreate(&start));
	HANDLE_ERROR(cudaEventCreate(&stop));
	HANDLE_ERROR(cudaEventRecord(start,0));

	//calculate force and send force components to be summed
	force_kernel<<<BlocksTet,Threads_Per_Block>>>(data->dev_dF
											  , data->dev_dFpitch
											  , data->dev_TetNodeRank
											  , Ntets 
											  , data->dev_v
											  , data->dev_vpitch
											  , data->dev_pe
											  , data->dev_TetVol
											  , data->dev_TetToNode
											  , data->dev_TetToNodepitch
                        , data->dev_swell
											  , dt*float(iKern));

	//sync threads before updating
	cudaThreadSynchronize();

	//sum forces and update positions
	updateKernel<<<BlocksNode,Threads_Per_Block>>>( data->dev_dF
												, data->dev_dFpitch
												, data->dev_F
												, data->dev_Fpitch
												, Nnodes
												, data->dev_nodeRank
												, data->dev_v
												, data->dev_vpitch
												, data->dev_r
												, data->dev_rpitch
												, data->dev_m);

	

	HANDLE_ERROR(cudaEventRecord(stop, 0));
	HANDLE_ERROR(cudaEventSynchronize(stop));
	HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));
	HANDLE_ERROR( cudaEventDestroy( start ));
	HANDLE_ERROR( cudaEventDestroy( stop ));
	
	

	//make sure kernel has finished before frame is printed 
	cudaThreadSynchronize(); 
	
	if((iKern)%iterPerFrame==0){
		
		//print calculation speed
		printf("\nIteration rate:  %f  iteartion/s \n kernel %d of %d\n"
								,1000.0/elapsedTime
								,iKern+1
								,NSTEPS);

	   //print frame
    printMetric(   data
						,host_data
						,Ntets
						,Nnodes
						,iKern+1);

		printVTKframe(   data
						,host_data
						,Ntets
						,Nnodes
						,iKern+1);
		printf("time = %f seconds\n", float(iKern)*dt);


		//print energy
		getEnergy(	 data
					,host_data
					,Ntets
					,Nnodes
					,pE
					,kE );
		fprintf(Eout,"%f %f %f %f\n",float(iKern)*dt,pE,kE,pE+kE);
		fflush(Eout);

    //graph->update(data,host_data);
    //graph->calcGauss();
    //graph->printFrame(iKern+1);
    //graph->dumpState(iKern+1);
    
	}//if((iKern+1)%iterPerFrame==0)


	

	
	//reset global mutex
	 HANDLE_ERROR( cudaMemset( g_mutex, 0, sizeof(int) ) );

	
	}//iKern

	fclose(Eout);

	




	//===================================================================




};

#endif

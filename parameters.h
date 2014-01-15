#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

//standard simulation parameters
#define NSTEPS			7000000					    //total number of iterations
#define dt              0.00000005					    //timestep [s]
#define iterPerFrame    10000                            //iterations per printed frame



//convert mesh length scale to cm
#define meshScale        0.1         //--[ cm / mesh unit]                      

//confineing prameters
#define CONFTOP 0.7        //confined top
#define CONFBOTTOM 0.2     //confined bottom
#define CONFFORCE  100000.0     //confining force if above CONFTOP or below CONFBOTTOM


// define lamda and mu
#define LAMBDA 56000000.0  
#define MU   2280000.0

//Elasticity constants (Lame' Coefficients)
//  -there is a factor of two off here from 
//  -the physical values, not sure if it
//  -is a half or two... I will find out
#define cxxxx			 29140000.0	  //--[ g / cm * s^2 ]			
#define cxxyy			 28000000.0	  //--[ g / cm * s^2 ]			
#define cxyxy			 570000.0	  //--[ g / cm * s^2 ]	


//Density of elastomer material
#define materialDensity  1.2   //--[ g / cm^3 ]  =  [  10^3 kg / m^3 ]


//scalar velocity dampening
//each velocity multiplied by this at each step
#define damp             0.999      //--[ unitless ]


//x and y dimensions of n profile
//input arrays 
#define inX 200
#define inY 200

//Threads per block to exicute
//100 seems pretty optimal on GTX275
//might be better with larger on 
//better card/ differnt card
#define TPB				100			



//maximum number of tetrahedra
//a Node can belone to
#define MaxNodeRank     65							



//constants declared on the stack for speed
#define PI				3.14159265359
#define dt2o2           (dt*dt)/2.0					    //for speed
#define dto2             dt/2.0						    //for speed

#endif //__PARAMETERS_H__

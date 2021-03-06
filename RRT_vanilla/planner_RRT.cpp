/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include "tree.hpp"
#include <iostream>
#include <fstream>
//matrices

#include <random>
#include <cmath> //cos
#include <math.h>
#include <algorithm> //fill
#include <fstream>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <string>
#include <vector>
#include <valarray> 
#include <list>
/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]


/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654
using namespace std;

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

static void planner(
		   double*	map,
		   int x_size,
 		   int y_size,
        double* armstart_anglesV_rad,
        double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	//no plan by default
  int goalReached = 0;
	*plan = NULL;
	*planlength = 0;
  int NoOfRuns = 0;
 //  // int** parent;
 //  // Tree tree;

 //  Tree tree(armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs);
 //  goalReached = tree.BuildRRT(map, x_size, y_size);
 //  int len = 0;


 //  if(goalReached)
 //  {
 //    vertex* temp = new vertex;
  
 //    temp = tree.goal;
 //    vector<vector <double>> pathVec;
    
 //    vector<double> currVec;

 //    // Pushing the backtracked tree in a vector
 //    while (temp!= tree.start)
 //    {
 //      currVec = temp->theta;
 //      pathVec.push_back(currVec);
 //      len = len+1;
 //      mexPrintf("length of tree: %d \n", len);
 //      temp = temp->parent;
 //    }
 //    // pushing in start
 //    if(temp!=NULL){
 //      currVec = temp->theta;
 //      pathVec.push_back(currVec);
 //      len = len+1;
 //    }


 //    *plan = (double**) malloc(len*sizeof(double*));

 //    for (int i = 0; i<len; i++){
 //        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
 //        for(int j = 0; j < numofDOFs; j++){
 //            // *(*(*(plan+i)+j)+0) = tree[tree.size()-i-1][j];
 //          (*plan)[i][j]=pathVec[len-i-1][j];
 //        }
 //    }
 //  }
 //  *planlength = len;
  while(!goalReached){
    NoOfRuns = NoOfRuns + 1;
    goalReached = mainRun(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
  }
  mexPrintf("\n No. of runs = %d", NoOfRuns);
  return;
}

//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    
    /* Check for proper number of arguments */    
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Three input arguments required."); 
    } else if (nlhs != 2) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 

    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);

    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))){
        	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);
    
    //call the planner
    double** plan = NULL;
    int planlength = 0;
    mexPrintf("mex working fine ");
    
    planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
    
    printf("planner returned plan of length=%d\n", planlength); 
    
    /* Create return values */
    if(planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
                plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    int* planlength_out = (int*) mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    
    return;
    
}






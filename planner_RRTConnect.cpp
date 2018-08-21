/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include "planner.h"
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

float goalProbability = 0.0;

typedef struct {
  int X1, Y1;
  int X2, Y2;
  int Increment;
  int UsingYIndex;
  int DeltaX, DeltaY;
  int DTerm;
  int IncrE, IncrNE;
  int XIndex, YIndex;
  int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
  params->UsingYIndex = 0;

  if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
    (params->UsingYIndex)++;

  if (params->UsingYIndex)
    {
      params->Y1=p1x;
      params->X1=p1y;
      params->Y2=p2x;
      params->X2=p2y;
    }
  else
    {
      params->X1=p1x;
      params->Y1=p1y;
      params->X2=p2x;
      params->Y2=p2y;
    }

   if ((p2x - p1x) * (p2y - p1y) < 0)
    {
      params->Flipped = 1;
      params->Y1 = -params->Y1;
      params->Y2 = -params->Y2;
    }
  else
    params->Flipped = 0;

  if (params->X2 > params->X1)
    params->Increment = 1;
  else
    params->Increment = -1;

  params->DeltaX=params->X2-params->X1;
  params->DeltaY=params->Y2-params->Y1;

  params->IncrE=2*params->DeltaY*params->Increment;
  params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
  params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

  params->XIndex = params->X1;
  params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
  if (params->UsingYIndex)
    {
      *y = params->XIndex;
      *x = params->YIndex;
      if (params->Flipped)
        *x = -*x;
    }
  else
    {
      *x = params->XIndex;
      *y = params->YIndex;
      if (params->Flipped)
        *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
  if (params->XIndex == params->X2)
    {
      return 0;
    }
  params->XIndex += params->Increment;
  if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
    params->DTerm += params->IncrE;
  else
    {
      params->DTerm += params->IncrNE;
      params->YIndex += params->Increment;
    }
  return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
		   int x_size,
 		   int y_size)

{
	bresenham_param_t params;
	int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
    
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
            return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
		   int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
    
 	//iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
    y1 = 0;
	for(i = 0; i < numofDOFs; i++)
	{
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
				return 0;
	}    
}
double smallerAngle(double differ)
{
  if (differ > PI)
    return (differ - 2*PI);
  else if (differ < -PI)
    return (2*PI + differ);
  else 
    return differ;
}

float distNorm(double* q, double* qrand1, int numofDOFs)
{
  // float p = 2;
  double temp1 = 0;
  double dist = 0;
  vector<double> diff(numofDOFs,0);
  for(int i =0;i<numofDOFs; i++){
    // diff[i] = (*q)[i] - (*qrand1)[i];
    temp1 = smallerAngle(*(q+i) - *(qrand1+i));
    diff[i] = temp1*temp1;
    dist = dist + diff[i];
  }

  dist = std::pow(dist, 1.0/(float) numofDOFs);
  // mexPrintf("Distby =%f  ", dist);
  return dist;
}


double* RandRRT(double* qgoal, float goalProbability, int numofDOFs, double*  map, int x_size, int y_size){
  // srand (static_cast <unsigned> (time(0)));
  double* qrand = (double*) malloc(numofDOFs*sizeof(double));
  //srand (time(NULL));
  unsigned startSeed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(startSeed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    //float r = distribution(generator);
  //float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    double d = distribution(generator);
  if (d<goalProbability){
    // cout<<"r: "<<r<<= %endl;
     // mexPrintf("\ngoalProbability= %f\n", d);
    qrand=qgoal;
  }
  else{
    for(int i =0;i<numofDOFs; i++){
      *(qrand+i) = 2*PI*(distribution(generator));
    }
    //qrand = randu<fcolvec>(qgoal.n_elem)*M_PI;
  }
  if (!IsValidArmConfiguration(qgoal, numofDOFs, map, x_size, y_size))
  {
    mexPrintf("\nRandom regeneration---------------------\n");
    qrand = RandRRT(qgoal, goalProbability, numofDOFs, map, x_size, y_size);
  }

  return qrand;
}

//Replace this for RRT*
int extendRRT(state *qnear, state* qrand, float step_size, int numofDOFs, double*  map, int x_size, int y_size){ //Extending the tree to qnear after checking for collisions
  int flag;
  state *qnew1;
  qnew1 = new state;
  qnew1->q = (double*) malloc(numofDOFs*sizeof(double)); 
  //cout<<"Dist_norm_Pre?\n";
  float distby = distNorm(qrand->q, qnear->q, numofDOFs);
  double move = min(distby, step_size);

  vector<double> qdelta(numofDOFs,0);

  if (move == distby){
    qrand->PREV = qnear;
    qnear->NEXT.push_back(qrand);
    flag = 2; //reached
    return flag;
  }
  else{
    for(int i =0;i<numofDOFs; i++){
      qdelta[i] = (*(qrand->q+i) - *(qnear->q+i))*(move/distby);
      *(qnew1->q+i) = *(qnear->q+i) + qdelta[i];
    }

    if (IsValidArmConfiguration(qnew1->q, numofDOFs, map, x_size, y_size))
    {
      // mexPrintf("\nExtending----------------------\n");
      qnew1->PREV = qnear;
      qnear->NEXT.push_back(qnew1);
      flag = 1; //Extended
      // for(int i = 0 ; i< numofDOFs; i++){
      //   mexPrintf("qnew_inside: %f \n", *(qnew1->q+i));
      // }
      return flag;
    }

    else
    {
      mexPrintf("\nCollision detected\n");
      flag = 0; //Stuck
      return flag;
      //cout<<"Collision detected";
      //qnew1 = NULL;
    }
  }

  // return flag;
}

state *connect(state *qnearGoal, state *qnew, float step_size, int numofDOFs, double*  map, int x_size, int y_size){
  
  state *preserveConnecion;
  preserveConnecion = new state;
  preserveConnecion->q = (double*) malloc(numofDOFs*sizeof(double));
  preserveConnecion = qnew->PREV;
  int flag = 1;
  while(flag == 1){
    flag = extendRRT(qnearGoal, qnew, step_size, numofDOFs, map, x_size, y_size);
    if(flag == 1){
      
      qnearGoal = qnearGoal->NEXT[(qnearGoal->NEXT.size()-1)];
        // for(int i = 0 ; i< numofDOFs; i++){
        //   mexPrintf("\nqnearGoal: %f ", *(qnearGoal->q+i));
        // }
    }
  }
  if(flag ==2){
    qnew->PREV = preserveConnecion;
  }
  return qnearGoal;
}

state *nearestNeighborRRT(state *qi, double* qrand, float nearDist, int numofDOFs){ //Finding in the tree the nearest neighbour to qrand

  // Tree search (recursive)
  state *qnear;
  double dist2 = distNorm(qi->q, qrand, numofDOFs);
  // mexPrintf("Dist =%f  ", dist2);
  
  if (dist2 < nearDist)
  {
    qnear = qi;
    nearDist = dist2;
  }

   // mexPrintf("qnear_next size: %d \n", qi->NEXT.size());
  
  //cout<<"qi"<<qi->q<<"\n";
  if(qi->NEXT.size() == 0){

    // for(int i = 0 ; i< numofDOFs; i++){
    //   mexPrintf("qnear_inside: %f \n", *(qnear->q+i));
    // }
    return qnear;
  }
  else
  {
    for(int i = 0 ; i < qi->NEXT.size(); i++)
    {

      // qi = qi->NEXT[i];
      
      qnear = nearestNeighborRRT(qi->NEXT[i], qrand, nearDist, numofDOFs);
      // qi = qi->PREV;
    }
    
  }

  return qnear;
}

void printTreeBackward(state *temp, state *qgoal, state *qinit, int numofDOFs){
  
  vector<double> currVec(numofDOFs,0);

  while (temp != qgoal && temp != NULL && temp != qinit)
  {
    for(int i = 0;i<numofDOFs; i++){
      currVec[i] = *(temp->q+i);
      mexPrintf("currVec = %f ",currVec[i]);
    }
    mexPrintf("\n");
    temp = temp->PREV;
  }
  if(temp!=NULL){
    for(int i = 0;i<numofDOFs; i++){
        currVec[i] = *(temp->q+i);
        mexPrintf("CurrVec = %f ",currVec[i]);
    }
  }
}

void swap(state *s1, state *s2){
  state *temp;
  temp = s1;
  s1 = s2;
  s2 = temp;
  return;
}
state *BuildRRTConnect(state *qinit, state *qgoal, int maxconfigs, float step_size, double*  map, int x_size, int y_size, int numofDOFs){
  int flag;
  state *qnear;
  qnear = new state;
  qnear->q = (double*) malloc(numofDOFs*sizeof(double));

  state *qnew;
  qnew = new state;
  qnew->q = (double*) malloc(numofDOFs*sizeof(double));

  state *qrand;
  qrand = new state;
  qrand->q = (double*) malloc(numofDOFs*sizeof(double));

  state *qtemp;
  qtemp = new state;
  qtemp->q = (double*) malloc(numofDOFs*sizeof(double));

  state *qnearGoal;
  qnearGoal = new state;
  qnearGoal->q = (double*) malloc(numofDOFs*sizeof(double));

  state *qconnect;
  qconnect = new state;
  qconnect->q = (double*) malloc(numofDOFs*sizeof(double));

  // mexPrintf("\nChecking qrand\n");
  int itr;
  for (itr = 1; itr <=100; itr++){
         
    float nearDist = 100000;
    mexPrintf("\nConfig no.= %d", itr);
    qrand->q = RandRRT(qgoal->q, goalProbability, numofDOFs, map, x_size, y_size);
    // for(int i = 0 ; i< numofDOFs; i++){
    //   mexPrintf("\nqrand: %f ", *(qrand->q+i));
    // }
    
    qnear = nearestNeighborRRT(qinit, qrand->q, nearDist, numofDOFs);
    // for(int i = 0 ; i< numofDOFs; i++){
    //   mexPrintf("\nqnear: %f ", *(qnear->q+i));
    // }
    flag = extendRRT(qnear, qrand, step_size, numofDOFs, map, x_size, y_size);

    mexPrintf("Flag %d \n ", flag);
    if(flag == 1 || flag == 2){
      nearDist = 100000;
      qnew = qnear->NEXT[(qnear->NEXT.size()-1)];
      // for(int i = 0 ; i< numofDOFs; i++){
      //   mexPrintf("\nqnew: %f ", *(qnew->q+i));
      // }
      qnearGoal = nearestNeighborRRT(qgoal, qnew->q, nearDist, numofDOFs);
      // for(int i = 0 ; i< numofDOFs; i++){
      //   mexPrintf("\nqnearGoal: %f ", *(qnearGoal->q+i));
      // }
      qconnect = connect(qnearGoal, qnew, step_size, numofDOFs, map, x_size, y_size);
      for(int i = 0 ; i< numofDOFs; i++){
        mexPrintf("\nqconnect: %f ", *(qconnect->q+i));
      }
      if (distNorm(qconnect->q, qnew->q, numofDOFs)< step_size)
      {
        mexPrintf("Goal reached at %d\n ", itr);
        printTreeBackward(qconnect, qgoal, qinit, numofDOFs);
        printTreeBackward(qnew, qgoal, qinit, numofDOFs);
        qnew->NEXT.push_back(qconnect);
        for(int i = 0 ; i< numofDOFs; i++){
          mexPrintf("\nqnew: %f ", *(qnew->q+i));
        }
        
        return qnew;
      }
      mexPrintf("backward \n");
      printTreeBackward(qconnect, qgoal, qinit, numofDOFs);
      mexPrintf("forward \n");
      printTreeBackward(qnew, qgoal, qinit, numofDOFs);
      qtemp = qinit;
      qinit = qgoal;
      qgoal = qtemp;
    }

  }


  return NULL;
}

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
	*plan = NULL;
	*planlength = 0;
  int** parent;
  
    //for now just do straight interpolation between start and goal checking for the validity of samples
    //but YOU  WILL WANT TO REPLACE THE CODE BELOW WITH YOUR PLANNER
  // double step_size = PI/20;
  double step_size = 0.01;
  int maxconfigs = 10000;
  state *qinit;
  qinit = new state;
  qinit->q = (double*) malloc(numofDOFs*sizeof(double));
  qinit->q = armstart_anglesV_rad;

  //Defining goal state
  state *qgoal;
  qgoal = new state;
  qgoal->PREV = NULL;
  qgoal->q = (double*) malloc(numofDOFs*sizeof(double));
  qgoal->q = armgoal_anglesV_rad;

  state *qmid;
  qmid = new state;
  qmid->q = (double*) malloc(numofDOFs*sizeof(double));
  
  
  // BuildRRT(qinit, qgoal, maxconfigs, step_size, map, x_size, y_size, numofDOFs);
  qmid = BuildRRTConnect(qinit, qgoal, maxconfigs, step_size, map, x_size, y_size, numofDOFs);
    for(int i = 0 ; i< numofDOFs; i++){
      mexPrintf("\nqmid: %f ", *(qmid->q+i));
    }
  // mexPrintf("\nSize: %d ", qmid->NEXT.size());
    for(int i = 0 ; i< numofDOFs; i++){
      mexPrintf("qstart: %f \n", *(qinit->q+i));
    }
    for(int i = 0 ; i< numofDOFs; i++){
      mexPrintf("qgoal: %f \n", *(qgoal->q+i));
    }
  vector <vector <double>> tree;
  vector <vector <double>> tree1;
  vector <vector <double>> tree2;
  state *temp;
  temp = new state;
 

  vector<double> currVec(numofDOFs,0);

//Reading from qmid to qinit
  if(qmid!=NULL){
    temp = qmid;
  }
  else{
    temp = NULL;
  }

  int len1 = 0;

  while (temp != qinit && temp != NULL && temp != qgoal)
  {
    for(int i = 0;i<numofDOFs; i++){
      currVec[i] = *(temp->q+i);
      mexPrintf("currVec = %f ",currVec[i]);
    }
    mexPrintf("\n");
    tree1.push_back(currVec);
    len1 = len1+1;
    
    temp = temp->PREV;
  }
  
  if(temp!=NULL){
    for(int i = 0;i<numofDOFs; i++){
        currVec[i] = *(temp->q+i);
        mexPrintf("CurrVec = %f ",currVec[i]);
    }
    tree1.push_back(currVec);
    len1 = len1+1;
  }
  
  // return;
  mexPrintf("\n\nQmid to qmid\n\n----------------------------------- ");
//Reading from qmid to qinit
  if(qmid!=NULL){
    if(qmid->NEXT.size()!=0){
      temp = qmid->NEXT[(qmid->NEXT.size()-1)];
              for(int i = 0 ; i< numofDOFs; i++){
          mexPrintf("\ntemp %f ", *(temp->q+i));
        }
    }
  }
  else{
    temp = NULL;
  }
 
  int len2 = 0;
  while (temp != qgoal && temp != NULL && temp != qinit)
  {
    for(int i = 0;i<numofDOFs; i++){
      currVec[i] = *(temp->q+i);
      mexPrintf("currVec = %f ",currVec[i]);
    }
    mexPrintf("\n");
    tree2.push_back(currVec);
    len2 = len2+1;
    temp = temp->PREV;
  }
  if(temp!=NULL){
    for(int i = 0;i<numofDOFs; i++){
        currVec[i] = *(temp->q+i);
        mexPrintf("CurrVec = %f ",currVec[i]);
    }
    tree2.push_back(currVec);
    len2 = len2+1;
  }

  *plan = (double**) malloc((len1+len2)*sizeof(double*));
  
  mexPrintf("\nlen1: %d ", len1);
  mexPrintf("\nlen2: %d ", len2);
 
if(len1>0 && len2>0){


  if(tree1[len1-1][0] == *(armstart_anglesV_rad) && tree1[len1-1][1] == *(armstart_anglesV_rad+1))
  {

    for(int i = 0; i<len1; i++){
      (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
      for(int j = 0; j < numofDOFs; j++){
        (*plan)[i][j]=tree1[len1-i-1][j];
      }
    }
    for(int i = 0; i<len2; i++){
      (*plan)[len1+i] = (double*) malloc(numofDOFs*sizeof(double));
      for(int j = 0; j < numofDOFs; j++){
        (*plan)[len1+i][j]=tree2[i][j];
      }
    }

  }
  else{
    for(int i = 0; i<len2; i++){
      (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
      for(int j = 0; j < numofDOFs; j++){
        (*plan)[i][j]=tree2[len2-i-1][j];
      }
    }
    for(int i = 0; i<len1; i++){
      (*plan)[len2+i] = (double*) malloc(numofDOFs*sizeof(double));
      for(int j = 0; j < numofDOFs; j++){
        (*plan)[len2+i][j]=tree1[i][j];
      }
    }
  }
}

*planlength = len1+len2;


/*    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf)
        {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = numofsamples;*/
    
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






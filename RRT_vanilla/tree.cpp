#include <math.h>
#include "mex.h"
#include "tree.hpp"
#include <ctime>
#include <chrono>
#include <random>

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654
using namespace std;

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

//Pre-defined functions
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


int IsValidLineSegment(double x0, double y0, double x1, double y1, double*  map,
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

int IsValidArmConfiguration(double* angles, int numofDOFs, double*  map,
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

void print_vector(vector<double> vecToPrint, int numofDOFs)
{
  for(int j = 0; j < numofDOFs; j++)
  {
    mexPrintf("%f    ",vecToPrint[j]);
  }
  mexPrintf("\n");
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

vector<double> Tree::RandRRT(double*  map, int x_size, int y_size){

  vector<double> qrand;
  unsigned startSeed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(startSeed);
  uniform_real_distribution<double> distribution(0.0,1.0);
  
  double d = distribution(generator);

  if (d<goalProbability){
    qrand = goal->theta;
  }
  else{
    for(int i =0;i<numofDOFs; i++){
      qrand.push_back(2*PI*(distribution(generator)));
    }
  }
  
  double *check = new double[numofDOFs];

  for(int i = 0; i < numofDOFs; i++)
  {
    *(check + i) = goal->theta[i];
  }

  if (!IsValidArmConfiguration(check, numofDOFs, map, x_size, y_size))
  {
    qrand = RandRRT(map, x_size, y_size);
  }

  return qrand;
}

void Tree::nearestNeighborRRT(vector<double> qrand, vertex **nearestVertex){ //Finding in the tree the nearest neighbour to qrand

  double shortestDistance = 10000.0;
  double distance = 0;
  int i = 0;
  double difference[numofDOFs];

  for (list<vertex*>::iterator it= vertices.begin(); it != vertices.end(); ++it)
  {
    distance = 0;

    for (i = 0; i < numofDOFs; ++i)
    {
      distance += pow(smallerAngle((qrand[i] - (*it)->theta[i])),2);
    }

    distance = pow(distance,0.5);
    
    if (distance < shortestDistance)
    {
      shortestDistance = distance;
      *nearestVertex = (*it);
    }
  }

}

int Tree::extendRRT(vector<double> qrand, double*  map, int x_size, int y_size){ //Extending the tree to qnear after checking for collisions

  int flag;
  double distBwVectors = 0;
  vector<double> diffVector;
  double distToGoal = 0;
  double nextAngle;
  vector<double> qnew;
  vertex *nearestVertex = new vertex;
  vertex *newVertex = new vertex;
  Tree::nearestNeighborRRT(qrand, &nearestVertex);
  double *check = new double[numofDOFs];
  // mexPrintf("\n");
  // mexPrintf("Random: ");  
  // print_vector(qrand, numofDOFs);
  // mexPrintf("Nearest: "); 
  // print_vector(nearestVertex->theta, numofDOFs);
  
  for (int i = 0; i < numofDOFs; ++i)
  {
    diffVector.push_back(smallerAngle(qrand[i] - nearestVertex->theta[i]));
    distBwVectors += pow(diffVector[i],2);
  }
  distBwVectors = pow(distBwVectors, 0.5);
  // mexPrintf("Dist to goal: %f\n", distToGoal);
  

  if(distBwVectors > step_size)
  {
    for (int i = 0; i < numofDOFs; ++i)
    {
      nextAngle = nearestVertex->theta[i] + step_size*(diffVector[i])/distBwVectors;
      // mexPrintf("Next angle test = %g  ", nextAngle);
      if (nextAngle < 0)
        qnew.push_back((2*PI + nextAngle));
      else if (nextAngle > 2*PI)
        qnew.push_back((nextAngle - 2*PI));
      else
        qnew.push_back((nextAngle));
    }
    
  }
  else
  {
    for (int i = 0; i < numofDOFs; ++i)
      qnew.push_back((qrand[i]));
  }

  for(int i = 0; i < numofDOFs; i++)
  {
    *(check + i) = qnew[i];
  }

  if (IsValidArmConfiguration(check, numofDOFs, map, x_size, y_size))
  {
    distToGoal = 0;
    //adding new vertex to the tree;
    // mexPrintf("adding new vertex to the tree\n");
    // mexPrintf("NewVertex: ");
    // print_vector(qnew, numofDOFs);
    newVertex->theta = qnew;
    newVertex->parent = nearestVertex;
    nearestVertex->children.push_back(newVertex);
    vertices.push_back(newVertex);

    for (int i = 0; i < numofDOFs; ++i)
    {
      distToGoal += pow(smallerAngle(newVertex->theta[i] - goal->theta[i]),2);
    }
    distToGoal = pow(distToGoal,0.5);
    if(distToGoal < step_size){
      flag = 1; //goal reached
      goal->parent = newVertex;
      newVertex->children.push_back(goal);
    }
    else{
      flag = 2; //Advanced
    }
    
  }
  else{
    flag = 3; //Blocked
  }
  return flag;
}

int Tree::BuildRRT(double* map, int x_size, int y_size){

  vector<double> qrand;

  int i;

  for (i = 1; i <=maxconfigs; i++){
         
    float nearDist = 100000;
    qrand = RandRRT(map, x_size, y_size);
    
    int flag = extendRRT(qrand, map, x_size, y_size);
    
    if(flag == 1){
      mexPrintf(" GOAL REACHED!! \n");
      return 1; //Goal reached
    }
  }
  mexPrintf(" GOAL NOT REACHED!! \n");
  return 0;
}

int mainRun(double*  map, int x_size, int y_size, double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs, double*** plan, int* planlength){
  *plan = NULL;
  *planlength = 0;
  int goalReached;
  Tree tree(armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs);
  goalReached = tree.BuildRRT(map, x_size, y_size);
  int len = 0;

  vertex* temp = new vertex;

  temp = tree.goal;
  vector<vector <double>> pathVec;
  
  if(goalReached)
  {
    vertex* temp = new vertex;
  
    temp = tree.goal;
    vector<vector <double>> pathVec;
    
    vector<double> currVec;

    // Pushing the backtracked tree in a vector
    while (temp!= tree.start)
    {
      currVec = temp->theta;
      pathVec.push_back(currVec);
      len = len+1;
      temp = temp->parent;
    }
    // pushing in start
    if(temp!=NULL){
      currVec = temp->theta;
      pathVec.push_back(currVec);
      len = len+1;
    }


    *plan = (double**) malloc(len*sizeof(double*));

    for (int i = 0; i<len; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double));
        for(int j = 0; j < numofDOFs; j++){
            // *(*(*(plan+i)+j)+0) = tree[tree.size()-i-1][j];
          (*plan)[i][j]=pathVec[len-i-1][j];
        }
    }
  }
  *planlength = len;
  return goalReached;
}
// void Tree::getPathFromTree(std::list<int> * pathVertices, int goalID)
// {
//   pathVertices->push_back(goalID);
//   int parentVertex = goalID; 
//   while (parentVertex!=1)
//   {
//     findParentVertex(pathVertices->back(), &parentVertex);
//     pathVertices->push_back(parentVertex);
//   }
// }

// void Tree::print_vertex(vertex *currentNode, int numofDOFs)
// {
//   int i = 0;
//   mexPrintf("vertex is: \n");
//   for (i = 0; i < numofDOFs; ++i)
//   {
//     mexPrintf("%g    ",currentNode->theta[i]);
//   }
//   mexPrintf("\n ID is: %d\n \n", currentNode->vertexID);
// }


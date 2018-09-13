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
#define MAX(A, B) ((A) > (B) ? (A) : (B))
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
// 2 (contact layer) and 0 (free space) are valid configurations; 1 (onstacle) and outside the map are invalid configurations
int IsValidState(double* state, int stateSize, double*  map, int x_size, int y_size){
  double x = state[0];
  double y = state[1];

  if(x < 0 || x >= x_size ||
  y < 0 || y >= y_size)
  return 0;
  short unsigned int nX, nY;
  ContXY2Cell(x, y, &nX, &nY, x_size, y_size);
  if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
        return 0;
  return 1;
}

int IsInContactLayer(double* state, int stateSize, double*  map,
  int x_size, int y_size){
  double x = state[0];
  double y = state[1];

  if(x < 0 || x >= x_size ||
  y < 0 || y >= y_size)
  return 0;
  short unsigned int nX, nY;
  ContXY2Cell(x, y, &nX, &nY, x_size, y_size);
  if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 2)
        return 1;
  return 0;
}

void print_vector(VectorXd vecToPrint, int stateSize)
{
  for(int j = 0; j < stateSize; j++)
  {
    mexPrintf("%f    ",vecToPrint[j]);
  }
  mexPrintf("\n");
}

void print_particleMatrix(MatrixXd matToPrint, int stateSize, int numofParticles)
{
  // mexPrintf("Can i print \n");
  for(int i = 0; i < stateSize; i++)
  {
    for (int j = 0; j< numofParticles; j++){
      mexPrintf("%f  ",matToPrint(i,j));
      // mexPrintf("Can i attempt to print \n");
    }
    mexPrintf("\n");
  }
  mexPrintf("\n");
}

double* toDoubleVector(VectorXd q, int stateSize){
  double *check = new double[stateSize];
  for(int i = 0; i < stateSize; i++)
  {
    *(check + i) = q(i);
  }
  return check;
}
//Convert using Armandillo
double Tree::costCalc(MatrixXd particleMatrix, VectorXd qrand)
{
  //-----------------------------------------------------------------------------------------------------------------------------------------------//
  //------Normalize both the distances to the interval [0,1] by dividing them by the maximum observed value over all the samples-------------------//
  //-----------------------------------------------------------------------------------------------------------------------------------------------//

  int p = 0;
  double totCost;
  double CovCost; //Cost resembling the standard devaition of the particle set
  double distCost; //Cost resembling the distance of rand from the mean of the particle set
  
  MatrixXd deltaMatrix(stateSize,numofParticles);
  VectorXd meanParticle(stateSize);
  MatrixXd covMat(stateSize,stateSize);

  meanParticle = particleMatrix.rowwise().mean();
  deltaMatrix = particleMatrix.colwise() - meanParticle;
  covMat = (deltaMatrix * deltaMatrix.transpose())/numofParticles;
  CovCost = sqrt(covMat.trace());


  //Finding the distance of rand from the mean of the particle set
  for (int i = 0; i < stateSize; i++){
    distCost += pow((qrand[i] - meanParticle(i)),2);
  }
  distCost = pow(distCost,0.5);
  totCost = gamma*CovCost + (1-gamma)*distCost;

  return totCost;
}


VectorXd Tree::RandRRT(double* map, int x_size, int y_size){

  VectorXd qrand(stateSize);
  unsigned startSeed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(startSeed);
  uniform_real_distribution<double> distribution(0.0,1.0);
  
  double d = distribution(generator);
  // mexPrintf("Distribution generator for qrand: %f\n", d);
  if (d<goalProbability){
    qrand = goal->particleMatrix.col(0);
  }
  else{
    //random value in the range of map

    qrand(0) = static_cast<double>(x_size)*(distribution(generator));
    qrand(1) = static_cast<double>(y_size)*(distribution(generator));
  }

  // mexPrintf("Printing typecast xsize: %f \n", static_cast<double>(x_size));
  if (!IsValidState(toDoubleVector(qrand, stateSize), stateSize, map, x_size, y_size))
  {
    qrand = RandRRT(map, x_size, y_size);
  }
  // mexPrintf("Printing qrand in function:\n");
	// print_vector(qrand, STATE_SIZE);
  return qrand;
}

void Tree::nearestNeighborRRT(VectorXd qrand, vertex **nearestVertex){ //Finding in the tree the nearest neighbour to qrand

  double shortestDistance = 10000.0;
  double maximumDistance = 0.0;
  double maximumCovariance = 0.0;
  double costTot = 0.0;

  VectorXd distVector;
  VectorXd covVector;
  distVector = VectorXd::Zero(vertices.size());
  covVector = VectorXd::Zero(vertices.size());

  double maxDist;
  double maxCov;
  //Normalized values
  double CovCost; //Cost resembling the standard devaition of the particle set
  double distCost = 0.0; //Cost resembling the distance of rand from the mean of the particle set
  
  MatrixXd deltaMatrix(stateSize,numofParticles);
  VectorXd meanParticle(stateSize);
  MatrixXd covMat(stateSize,stateSize);
  int temp = 0;
  for (list<vertex*>::iterator it1= vertices.begin(); it1 != vertices.end(); ++it1)
  {
  	
  	distCost = 0.0;
  	meanParticle = ((*it1)->particleMatrix).rowwise().mean();
  	deltaMatrix = ((*it1)->particleMatrix).colwise() - meanParticle;
	covMat = (deltaMatrix * deltaMatrix.transpose())/numofParticles;
	covVector[temp] = sqrt(covMat.trace());

	//Finding the distance of rand from the mean of the particle set
	for (int i = 0; i < stateSize; i++){
	  distCost += pow((qrand[i] - meanParticle[i]),2);
	}
	distVector(temp) = pow(distCost,0.5);

    if (distVector[temp] > maximumDistance)
    {
      maximumDistance = distVector[temp];
    }

    if (covVector[temp] > maximumCovariance)
    {
      maximumCovariance = covVector[temp];
    }
    temp = temp +1;
  }

  // mexPrintf("distVector filled: \n");
  // print_vector(distVector, vertices.size());

  // mexPrintf("covVector filled: \n");
  // print_vector(covVector, vertices.size());
  // mexPrintf("Printing max dist:%f\n", maximumDistance);
  // mexPrintf("Printing max cov:%f\n", maximumCovariance);
  temp = 0;
  for (list<vertex*>::iterator it2= vertices.begin(); it2 != vertices.end(); ++it2)
  {
    
    // mexPrintf("Printing vertices:\n");
    // print_particleMatrix((*it2)->particleMatrix, STATE_SIZE, NUMBEROFPARTICLES);
    costTot = gamma*covVector[temp]/maximumCovariance + (1-gamma)*distVector[temp]/maximumDistance;
    // mexPrintf("Printing tot dist:%f\n", costTot);
    if (costTot < shortestDistance)
    {
      shortestDistance = costTot;
      (*nearestVertex) = (*it2);
    }
    temp = temp +1;
  }

  return;
}

// void Tree::motion_model(){
//   return;
// }

int Tree::selectInput(vertex *nearestVertex){
  int action;
  unsigned startSeed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(startSeed);
  uniform_real_distribution<double> distribution(0.0,1.0);
  double g = distribution(generator);

  if(nearestVertex->inContact == 0 && g < gamma){
    action = 2; //Guarded: Come in contact along the direction of target
  }
  else if(nearestVertex->inContact == 0 && g > gamma){
    action = 1; //Connect
  }
  else if(nearestVertex->inContact == 1 && g < gamma){
    action = 3; //Slide
  }
  else if(nearestVertex->inContact == 1 && g > gamma){
    action = 1; //Connect
  }
  return action;
}

 
//moveToTargetOrContact: Returns 1 if contact occurs, returns 0 if target(qrand) is reached
//qrand is a Valid Configuration
int Tree::moveToTargetOrContact(VectorXd qrand, VectorXd *qnew, double*  map, int x_size, int y_size){
  int advance = 1;
  double cell_size = 1.0;
  double step = step_size;
  VectorXd q(stateSize);
  // int count = 0;

  while(advance == 1){
    if((qrand - *qnew).norm() < step_size){
      *qnew = qrand;
      advance = 0;
      // mexPrintf("TARGET REACHED\n");
      return advance;
    }
    else{
      q = *qnew + step*(qrand - *qnew)/(qrand - *qnew).norm(); //step towards qrand
      if(IsValidState(toDoubleVector(q, stateSize), stateSize, map, x_size, y_size)){
        // mexPrintf("STEPPED to location: ");
        // print_vector(q, STATE_SIZE);
        short unsigned int nX, nY;
        ContXY2Cell(q(0), q(1), &nX, &nY, x_size, y_size);
        // mexPrintf("STEPPED to cell= %d, %d\n", nX, nY);
        *qnew = q;
      }
      else{
        step = step - cell_size;
        // mexPrintf("Step reduced to=%f\n", step);
        if(step < cell_size){
          advance = 1;
          // mexPrintf("CONTACT MADE in move to contact function as step=%f\n", step);
          short unsigned int nX, nY;
          ContXY2Cell((*qnew)(0), (*qnew)(1), &nX, &nY, x_size, y_size);
          // mexPrintf("We are presumably in contact layer at %d, %d and is it in contact layer? = %d\n", nX, nY, IsInContactLayer(toDoubleVector(*qnew, stateSize), stateSize, map, x_size, y_size));
          return advance;
        }
      }
    }
    // count = count +1;
    // if(count == 20){
    //   mexPrintf("Count limit reached\n");
    //   break;
    // }
  }
  return 0;
}

// moveToContact: Connects to the nearest contact (obstacle/wall)
void Tree::moveToContact(VectorXd qrand, VectorXd *qnew, double*  map, int x_size, int y_size){
  int advance = 1;
  // int count = 0;
  double cell_size = 1.0;
  VectorXd qnew_init = *qnew;
  double step = step_size;
  VectorXd q(stateSize);
  while(advance == 1){
    q = *qnew + step*(qrand - qnew_init)/(qrand - qnew_init).norm(); //step towards qrand
    if(IsValidState(toDoubleVector(q, stateSize), stateSize, map, x_size, y_size)){
      // mexPrintf("STEPPED to location: ");
      // print_vector(q, STATE_SIZE);
      short unsigned int nX, nY;
      ContXY2Cell(q(0), q(1), &nX, &nY, x_size, y_size);
      // mexPrintf("STEPPED to cell= %d, %d\n", nX, nY);
      *qnew = q;
    }
    else{
      step = step - cell_size;
      // mexPrintf("Step reduced to=%f\n", step);
      if(step < cell_size){
        // mexPrintf("CONTACT MADE in move to contact function as step=%f\n", step);
        short unsigned int nX, nY;
        ContXY2Cell((*qnew)(0), (*qnew)(1), &nX, &nY, x_size, y_size);
        // mexPrintf("We are presumably in contact layer at %d, %d and is it in contact layer? = %d\n", nX, nY, IsInContactLayer(toDoubleVector(*qnew, stateSize), stateSize, map, x_size, y_size));
        advance = 0;
        return;
      }
    }
    // count = count +1;
    // if(count == 50){
    //   mexPrintf("Count limit reached\n");
    //   break;
    // }
  }
  return;
}

// Connect function: Move towards qrand till qrand is reached or contact occurs
// Returns 1 if it connects to qrand and returns 2 if it comes in contact with a wall/obstacle
void Tree::connect(MatrixXd targetMatrix, vertex *nearestVertex, vertex **newVertex, double*  map, int x_size, int y_size){
  // mexPrintf("Inside connect function:\n");
  int flag1 = 0; // is 1 if contact has occured in contact/target case
  int flag2 = 0;
  int temp1 = 0;
  int temp2 = 0;
  VectorXd qrand(stateSize);
  VectorXd qnew(stateSize);
  int setValidity = 1;
  (*newVertex)->inContact = 0;

  // Check particle set for validity
  for(int i = 0; i < numofParticles; i++){
    if (!IsValidState(toDoubleVector(targetMatrix.col(i), stateSize), stateSize, map, x_size, y_size)){
      setValidity = 0;
    }
  }

  if (setValidity == 1){
    // mexPrintf("Qrand is Valid particle set\n");
    for(int i = 0; i < numofParticles; i++){
      // mexPrintf("Particle: %d\n", i);
      qrand = targetMatrix.col(i);
      qnew = nearestVertex->particleMatrix.col(i);
      temp1 = moveToTargetOrContact(qrand, &qnew, map, x_size, y_size); //Returns 1 if contact occurs
      flag1 = MAX(temp1, flag1);
      (*newVertex)->particleMatrix.col(i) = qnew;
    }
  }
  else{
    // mexPrintf("Qrand is Invalid particle set\n");
    flag2 = 1;
    for(int i = 0; i < numofParticles; i++){
      // mexPrintf("Particle: %d\n", i);
      qrand = targetMatrix.col(i);
      qnew = nearestVertex->particleMatrix.col(i);
      moveToContact(qrand, &qnew, map, x_size, y_size);
      (*newVertex)->particleMatrix.col(i) = qnew;
    }

  }

  //Setting the inContact parameter
  if(MAX(flag1, flag2)){
    (*newVertex)->inContact = 1;
  }

  (*newVertex)->parent = nearestVertex;
  nearestVertex->children.push_back((*newVertex));
  
  return;
}

void Tree::guarded(MatrixXd targetMatrix, vertex *nearestVertex, vertex **newVertex, double*  map, int x_size, int y_size){

  VectorXd qrand(stateSize);
  VectorXd qnew(stateSize);
  
  for(int i = 0; i < numofParticles; i++){
    qrand = targetMatrix.col(i);
    qnew = nearestVertex->particleMatrix.col(i);
    moveToContact(qrand, &qnew, map, x_size, y_size);
    (*newVertex)->particleMatrix.col(i) = qnew;
  }

  (*newVertex)->inContact = 1;

  (*newVertex)->parent = nearestVertex;
  nearestVertex->children.push_back((*newVertex));
  
  return;
}

//moveToTargetOrContact: Returns 1 if new contact occurs, returns 2 if contact is lost(wont occur for obstacles), returns 3 if target(qrand) is reached
//qrand is a Valid Configuration
void Tree::slide(MatrixXd targetMatrix, vertex *nearestVertex, vertex **newVertex, double*  map, int x_size, int y_size){
  double cell_size = 1.0;
  // mexPrintf("Inside slide function:\n");
  int flag1 = 0; // is 1 if contact has occured in contact/target case
  int flag2 = 0;
  int temp1 = 0;
  int temp2 = 0;
  VectorXd qrand(stateSize);
  VectorXd qnew(stateSize);
  int setValidity = 1;
  (*newVertex)->inContact = 0;

  // VectorXd q(stateSize);
  MatrixXd projectedTargetMatrix(stateSize, numofParticles);
  MatrixXd dotProduct(1, numofParticles);
  VectorXd normal(stateSize);

  // Project target matrix to the contact surface of the nearest vertex
  normal << MAP_normalx(nearestVertex->particleMatrix(0,0) ,nearestVertex->particleMatrix(1,0)), MAP_normaly(nearestVertex->particleMatrix(0,0) ,nearestVertex->particleMatrix(1,0));
  dotProduct = ((targetMatrix - nearestVertex->particleMatrix).cwiseProduct(normal.replicate<1,NUMBEROFPARTICLES>())).colwise().sum();
  projectedTargetMatrix = targetMatrix - (dotProduct.replicate<STATE_SIZE,1>()).cwiseProduct(normal.replicate<1,NUMBEROFPARTICLES>());
  
  // int advance = 3;
  // int targetReached = 0;
  // int newContact = 0;
  // int lostContact = 0;
  // pair<double,double> newContactSurface;
  // VectorXd pointOnNewContactSurf(stateSize);
  // VectorXd lostContactEdge(stateSize);

  // Check particle set for validity
  for(int i = 0; i < numofParticles; i++){
    if (!IsValidState(toDoubleVector(projectedTargetMatrix.col(i), stateSize), stateSize, map, x_size, y_size)){
      setValidity = 0;
    }
  }

  if (setValidity == 1){ //If it is a valid particle set directly move to it unless you get a contact in between
    // mexPrintf("Qrand is Valid particle set\n");
    for(int i = 0; i < numofParticles; i++){
      // mexPrintf("Particle: %d\n", i);
      qrand = projectedTargetMatrix.col(i);
      qnew = nearestVertex->particleMatrix.col(i);
      temp1 = moveToTargetOrContact(qrand, &qnew, map, x_size, y_size); //Returns 1 if contact occurs
      flag1 = MAX(temp1, flag1);
      (*newVertex)->particleMatrix.col(i) = qnew;
    }
  }
  else{ //If it is an invalid particle set move to the contact in the direction of target matrix
    // mexPrintf("Qrand is Invalid particle set\n");
    flag2 = 1;
    for(int i = 0; i < numofParticles; i++){
      // mexPrintf("Particle: %d\n", i);
      qrand = projectedTargetMatrix.col(i);
      qnew = nearestVertex->particleMatrix.col(i);
      moveToContact(qrand, &qnew, map, x_size, y_size);
      (*newVertex)->particleMatrix.col(i) = qnew;
    }

  }
  
  // for(int i = 0; i < numofParticles; i++){
  //   qrand = projectedTargetMatrix.col(i);
  //   qnew = nearestVertex->particleMatrix.col(i);
  //   double step = step_size;
    // int temp1 = 0;
    // int temp2 = 0;
    // int temp3 = 0;

    // while((temp1+temp2+temp3) == 0){ 
    //   if((qrand - qnew).norm() < cell_size){
    //     step = step - cell_size;
    //     qnew = qrand;
    //     temp1 = 1;
    //   }
    //   else{
    //     q = qnew + step*(qrand - qnew)/(qrand - qnew).norm(); //step towards qrand
    //     if(IsValidState(toDoubleVector(q, stateSize), stateSize, map, x_size, y_size) && IsInContactLayer(toDoubleVector(q, stateSize), stateSize, map, x_size, y_size) == 1){
    //       qnew = q;
    //     }
    //     else {
    //       step = step - cell_size;
    //       if(step < cell_size){
    //         if(!IsValidState(toDoubleVector(q, stateSize), stateSize, map, x_size, y_size)){
    //           temp2 = 1; // New contact
    //           // newContactSurface = pair(MAP_normalx(qnew(0), qnew(1)), MAP_normaly(qnew(0), qnew(1)));
    //           pointOnNewContactSurf = qnew;
    //         }
    //         else{
    //           lostContactEdge = qnew;
    //           temp3 = 1; //Lost contact
    //         }
    //       }
    //     }
    //   }
    // }

      //Setting the inContact parameter
	if(MAX(flag1, flag2)){
	  (*newVertex)->inContact = 1;
	}
	(*newVertex)->parent = nearestVertex;
	nearestVertex->children.push_back((*newVertex));
  	

  // if(newContact > 0){
  //   advance = 1;
  //   normal[0] = MAP_normalx(pointOnNewContactSurf[0], pointOnNewContactSurf[1]);
  //   normal[1] = MAP_normaly(pointOnNewContactSurf[0], pointOnNewContactSurf[1]);
  //   dotProduct = (((*newVertex)->particleMatrix - pointOnNewContactSurf.replicate<1,NUMBEROFPARTICLES>()).cwiseProduct(normal.replicate<1,NUMBEROFPARTICLES>())).colwise().sum();
  //   (*newVertex)->particleMatrix = (*newVertex)->particleMatrix - (dotProduct.replicate<STATE_SIZE,1>()).cwiseProduct(normal.replicate<1,NUMBEROFPARTICLES>());
  //   (*newVertex)->inContact = 1;
  // }
  // else if(lostContact > 0){
  //   advance = 2;
  //   // ----------------------- Include this later, how to do it in higher dimensions?-----------------------------
  //   // projectToLostContact(newVertex, newContactSurface);
  //   (*newVertex)->inContact = 0;
  // }

  // motion_model(newVertex);
  // project(newVertex);
  return;
}

//return 1 if goal is reached 
int Tree::extendRRT(VectorXd qrand, double*  map, int x_size, int y_size){ //Extending the tree to qnear after checking for collisions

  int flag = 0;

  vertex *nearestVertex = new vertex;
  vertex *newVertex = new vertex;
  
  newVertex->particleMatrix = MatrixXd::Zero(stateSize, numofParticles);

  // Finding the nearest vertex to qrand
  nearestNeighborRRT(qrand, &nearestVertex);
  // mexPrintf("Nearest neighbour contact state = %d\n", nearestVertex->inContact);
  // mexPrintf("Printing nearestVertex matrix: \n");
  // print_particleMatrix(nearestVertex->particleMatrix, STATE_SIZE, NUMBEROFPARTICLES);

  
  //Finding the action to take based on nearest node
  int action = selectInput(nearestVertex);

  mexPrintf("Action taken = %d\n", action);
  MatrixXd targetMatrix(stateSize, numofParticles);

  //Find qtarget = qrand + (qnearest - qnearestmean) to extend nearest vertex to:
  targetMatrix = qrand.replicate<1,NUMBEROFPARTICLES>() + (nearestVertex->particleMatrix - (nearestVertex->particleMatrix.rowwise().mean()).replicate<1,NUMBEROFPARTICLES>());
  // mexPrintf("Printing target matrix:\n");
  // print_particleMatrix(targetMatrix, STATE_SIZE, NUMBEROFPARTICLES);

  //Calling functions based on actions selected
  switch(action) {
    case 1 : connect(targetMatrix, nearestVertex, &newVertex, map, x_size, y_size); break;
    case 2 : guarded(targetMatrix, nearestVertex, &newVertex, map, x_size, y_size); break;
    case 3 : slide(targetMatrix, nearestVertex, &newVertex, map, x_size, y_size);
  }
  //newVertex added is a valid config so is added to the tree

  //Check distance to goal, return 1 if goal is reached
  // mexPrintf("Printing newVertex matrix:\n");
  // print_particleMatrix(newVertex->particleMatrix, STATE_SIZE, NUMBEROFPARTICLES);
  vertices.push_back(newVertex);
  if((newVertex->particleMatrix.rowwise().mean() - goal->particleMatrix.col(0)).norm() < step_size){
    flag = 1; //Goal reached
    goal->parent = newVertex;
    newVertex->children.push_back(goal);
    vertices.push_back(goal);
  }

  return flag;
}
//return 1 if goal is reached
int Tree::BuildRRT(double* map, int x_size, int y_size){

  VectorXd qrand;

  for (int i = 1; i <=maxconfigs; i++){
    mexPrintf(" EXTEND ITERATION: %d \n", i);
    
    qrand = RandRRT(map, x_size, y_size);
    
    // mexPrintf("Printing qrand:\n");
    // print_vector(qrand, STATE_SIZE);

    int flag = extendRRT(qrand, map, x_size, y_size);
    // mexPrintf(" Size of vertices list after extend=%d \n\n", vertices.size());

    if(flag == 1){
      mexPrintf(" GOAL REACHED!! \n");
      return 1; //Goal reached
    }
  }
  mexPrintf(" GOAL NOT REACHED!! \n");
  return 0;
}

int mainRun(double*  map, int x_size, int y_size, double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int stateSize, double*** plan, int* planlength){
  *plan = NULL;
  *planlength = 0;
  int goalReached;
  
  Tree tree(armstart_anglesV_rad, armgoal_anglesV_rad, stateSize, map, x_size, y_size);

  // //Testing validity of nodes
  // VectorXd testValidity(STATE_SIZE);
  // testValidity << 19, 6;
  // mexPrintf("Is start a valid state: %d\n", IsValidState(toDoubleVector(tree.start->particleMatrix.col(0), stateSize), stateSize, map, x_size, y_size));
  // mexPrintf("Is goal a valid state: %d\n", IsValidState(toDoubleVector(tree.goal->particleMatrix.col(0), stateSize), stateSize, map, x_size, y_size));  
  // mexPrintf("Is testValidity=%f %f a valid state: %d\n", testValidity(0), testValidity(1), IsValidState(toDoubleVector(testValidity, stateSize), stateSize, map, x_size, y_size));  
  

  //Finding map values
  short unsigned int nX, nY;
  double x = 23.0;
  double y = 5.0;
  ContXY2Cell(x, y, &nX, &nY, x_size, y_size);
  // mexPrintf("\nMap value at %f, %f = %f\n", x, y, map[GETMAPINDEX(nX,nY,x_size,y_size)]);
  
  // mexPrintf("Printing map:\n");
  // for(int i = 0; i < x_size; i++){
  //   for(int j = 0; j<y_size; j++){
  //     mexPrintf("%f  ", map[GETMAPINDEX(i,j,x_size,y_size)]);
  //   }
  //  mexPrintf("\n"); 
  // }
  
  goalReached = tree.BuildRRT(map, x_size, y_size);
  int len = 0;
  
  vertex* temp = new vertex;
  temp = tree.goal;
  vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd> > pathVec;
  
  if(goalReached)
  {
    vertex* temp = new vertex;
    temp = tree.goal;
    MatrixXd currParticleMatrix(STATE_SIZE,NUMBEROFPARTICLES);

    // Pushing the backtracked tree in a vector
    while (temp!= tree.start)
    {
      currParticleMatrix = temp->particleMatrix;
      pathVec.push_back(currParticleMatrix);
      len = len+1;
      temp = temp->parent;
    }

    // pushing in start
    if(temp!=NULL){
      currParticleMatrix = temp->particleMatrix;
      pathVec.push_back(currParticleMatrix);
      len = len+1;
    }


    *plan = (double**) malloc(len*sizeof(double*));
    for (int i = 0; i<len; i++){
      (*plan)[i] = (double*) malloc((NUMBEROFPARTICLES*STATE_SIZE)*sizeof(double));
      for(int j = 0; j < (NUMBEROFPARTICLES*STATE_SIZE); j++){
          // *(*(*(plan+i)+j)+0) = tree[tree.size()-i-1][j];
        (*plan)[i][j]=pathVec[len-i-1](j%STATE_SIZE, (int)(j/STATE_SIZE));
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
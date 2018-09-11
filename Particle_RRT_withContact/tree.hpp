/*
 * tree.hpp
 * Contains a class that defines tree structure
 * Tree node must contain a particle set to define the state and a contact pair
 * We would currently work with a single contact so one contact pair per state
 */

#ifndef TREE__
#define TREE__

#include <list>
#include <vector>
#include <ctime>
#include <chrono>
#include <random>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
#define PI 3.141592654
#define NUMBEROFPARTICLES 5

#define STATE_SIZE 2
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;

typedef Matrix<double, STATE_SIZE, NUMBEROFPARTICLES> MatrixPd;

//Removed this as using eigen matrices made more sense for matrix calculations
// typedef struct particle
// {
//     vector<double> xy; // for a 2D world and a point robot this would be (x,y)
//     // float weight; // Weight of this particle, maybe used later for resampling (not sure)
// } Particle;

typedef struct vertex{
    // list<Particle> ParticleSet; //Particle set for state of the robot
    // MatrixXd particleMatrix;
    MatrixXd particleMatrix;
    bool inContact; // Boolean to know if the robot is in contact or not, default not in contact
    std::list<vertex*> children;
    vertex* parent;
} vertex;

class Tree
{

public:
	std::list<vertex*> vertices;
	double step_size;
	int stateSize; //size of xy
	int numofParticles;
	double gamma;
	vertex* start = new vertex;
	vertex* goal = new vertex;
	float goalProbability;
	int maxconfigs;
  MatrixXd MAP_normalx;
  MatrixXd MAP_normaly;

	Tree(double* pointBot_start_state_xy,
	    double* pointBot_goal_state_xy,
	    int nDOFs, double* map, int x_size, int y_size){

    MAP_normalx = MatrixXd::Zero(x_size, y_size);
    MAP_normaly = MatrixXd::Zero(x_size, y_size);
		//Initializing start and goal
		start->parent = NULL;
		goal->parent = NULL;
		step_size = 4.0;
  	stateSize = nDOFs;
  	numofParticles = NUMBEROFPARTICLES;
  	goalProbability = 0.9;
  	maxconfigs = 10;
  	gamma = 0.5;

  	unsigned startSeed = std::chrono::system_clock::now().time_since_epoch().count();
  	default_random_engine generator(startSeed);
  	uniform_real_distribution<double> distribution(0.0,7.0);
    start->particleMatrix = MatrixXd::Zero(stateSize, numofParticles);
    start->inContact = 0;
    goal->particleMatrix = MatrixXd::Zero(stateSize, numofParticles);
    goal->inContact = 0;
  	for (int i = 0; i < stateSize; i++)
  	{
  		for(int j = 0 ; j < numofParticles; j++){
        start->particleMatrix(i,j) = distribution(generator)+ *(pointBot_start_state_xy+i);
        goal->particleMatrix(i,j) = *(pointBot_goal_state_xy+i);
      }
  	}
    
  	start->inContact = 0;
  	vertices.push_back(start);

    // Initializing map normal vectors at borders of obstacles
    //Verticle obstacle
    for(int i = 0; i <= 13; i++)
    {
      MAP_normalx(i,45) = -1.0;
      // MAP_normalx(i,18) = 1.0;
    }
    MAP_normaly(13,46) = -1.0;
    MAP_normaly(13,47) = -1.0;
    MAP_normaly(13,48) = -1.0;
    MAP_normaly(13,49) = -1.0;


    //Horizontal obstacle
    for(int i = 0; i <= 38; i++)
    {
      MAP_normaly(17,i) = 1.0;
      MAP_normaly(22,i) = -1.0;
    }
    MAP_normalx(18,38) = 1.0;
    MAP_normalx(19,38) = 1.0;
    MAP_normalx(20,38) = 1.0;
    MAP_normalx(21,38) = 1.0;

    // Initializing map normal vectors at walls

    for(int i = 0; i <= 17; i++) //left edge
    {
      MAP_normalx(i,0) = 1.0;
    }
    for(int i = 22; i <= 49; i++) //left edge
    {
      MAP_normalx(i,0) = 1.0;
    }
    for(int i = 0; i <= 45; i++) //top edge
    {
      MAP_normaly(0,i) = -1.0;
    }
    // for(int i = 18; i <= 49; i++) //top edge
    // {
    //   MAP_normaly(0,i) = -1.0;
    // }
    for(int i = 13; i <= 49; i++) //right edge
    {
      MAP_normalx(i,49) = -1.0;
    }
    for(int i = 0; i <= 49; i++) //bottom edge
    {
      MAP_normaly(49,i) = 1.0;
    }

  }
public:

  double costCalc(MatrixXd particleMatrix, VectorXd qrand);
	VectorXd RandRRT(double* map, int x_size, int y_size);
	void nearestNeighborRRT(VectorXd qrand, vertex **nearestVertex);
  int selectInput(vertex *nearestVertex);
  int moveToTargetOrContact(VectorXd qrand, VectorXd *qnew, double*  map, int x_size, int y_size);
  void  moveToContact(VectorXd qrand, VectorXd *qnew, double*  map, int x_size, int y_size);
  void  connect(MatrixXd targetMatrix, vertex *nearestVertex, vertex **newVertex, double*  map, int x_size, int y_size);
  void  guarded(MatrixXd targetMatrix, vertex *nearestVertex, vertex **newVertex, double*  map, int x_size, int y_size);
  void  slide(MatrixXd targetMatrix, vertex *nearestVertex, vertex **newVertex, double*  map, int x_size, int y_size);
	int 	extendRRT(VectorXd qrand, double*  map, int x_size, int y_size);
	int 	BuildRRT(double*  map, int x_size, int y_size);


};

//Pre-defined functions
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

void 	ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size);
void 	get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params);
void 	get_current_point(bresenham_param_t *params, int *x, int *y);
int 	get_next_point(bresenham_param_t *params);
int 	IsValidLineSegment(double x0, double y0, double x1, double y1, double*  map, int x_size, int y_size);
int 	IsValidArmConfiguration(double* angles, int numofDOFs, double*  map, int x_size, int y_size);
int 	mainRun(double*  map, int x_size, int y_size, double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs, double*** plan, int* planlength);
int IsValidState(double* state, int stateSize, double*  map, int x_size, int y_size);
int IsInContactLayer(double* state, int stateSize, double*  map, int x_size, int y_size);
//Helper functions

float distNorm(double* q, double* qrand1, int numofDOFs);
void 	print_vector(VectorXd vecToPrint, int numofDOFs);
void  print_particleMatrix(MatrixXd matToPrint, int stateSize, int numofParticles);
double* toDoubleVector(VectorXd q, int stateSize);


#endif // TREE__
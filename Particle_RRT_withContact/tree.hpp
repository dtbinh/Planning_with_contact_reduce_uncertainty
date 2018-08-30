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
#include <eigen3/Eigen/Dense>

using namespace Eigen;
#define PI 3.141592654
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;

//Removed this as using eigen matrices made more sense for matrix calculations
// typedef struct particle
// {
//     vector<double> xy; // for a 2D world and a point robot this would be (x,y)
//     // float weight; // Weight of this particle, maybe used later for resampling (not sure)
// } Particle;

typedef struct vertex{
    // list<Particle> ParticleSet; //Particle set for state of the robot
    MatrixXd particleMatrix(stateSize,numofParticles);
    bool inContact = 0; // Boolean to know if the robot is in contact or not, default not in contact
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
	float gamma;
	vertex* start = new vertex;
	vertex* goal = new vertex;
	float goalProbability;
	int maxconfigs;
  vector<vector<pair<double,double>>> MAP_normal;
	Tree(double* pointBot_start_state_xy,
	    double* pointBot_goal_state_xy,
	    int nDOFs, double* map, int x_size, int y_size){

		//Initializing start and goal
		start->parent = NULL;
		// start->vertexID = 1;
		goal->parent = NULL;
		step_size = 0.3;// PI*5/180; //% degrees increment
  	stateSize = nDOFs;
  	numofParticles = 20;
  	goalProbability = 0.90;
  	maxconfigs = 30000;
  	gamma = 0.5;

  	unsigned startSeed = std::chrono::system_clock::now().time_since_epoch().count();
  	default_random_engine generator(startSeed);
  	uniform_real_distribution<double> distribution(0.0,9.0);

  	for (int i = 0; i < numofParticles; i++)
  	{
  		for(int j = 0 ; j < stateSize; j++){
        start->particleMatrix(j,i) = distribution(generator)+ *(armstart_anglesV_rad+i);
        goal->particleMatrix(j,i) = *(armgoal_anglesV_rad+i);
      }
  	}

  	start->inContact = FALSE;
  	vertices.push_back(start);

    // Initializing map normal vectors at borders of obstacles

    for(int i = 0; i <= 12; i++)
    {
      MAP_normal[i][14] = (-1,0);
      MAP_normal[i][18] = (1,0);
    }
    MAP_normal[13][15] = (0,-1);
    MAP_normal[13][16] = (0,-1);
    MAP_normal[13][17] = (0,-1);
    for(int i = 0; i <= 8; i++)
    {
      MAP_normal[17][i] = (0,1);
      MAP_normal[22][i] = (0,-1);
    }
    MAP_normal[18][9] = (1,0);
    MAP_normal[19][9] = (1,0);
    MAP_normal[20][9] = (1,0);
    MAP_normal[21][9] = (1,0);

    // Initializing map normal vectors at walls

    for(int i = 0; i <= 17; i++) //left edge
    {
      MAP_normal[i][0] = (1,0);
    }
    for(int i = 22; i <= 49; i++) //left edge
    {
      MAP_normal[i][0] = (1,0);
    }
    for(int i = 0; i <= 14; i++) //top edge
    {
      MAP_normal[0][i] = (0,-1);
    }
    for(int i = 18; i <= 49; i++) //top edge
    {
      MAP_normal[0][i] = (0,-1);
    }
    for(int i = 0; i <= 49; i++) //right edge
    {
      MAP_normal[i][49] = (-1,0);
    }
    for(int i = 0; i <= 49; i++) //bottom edge
    {
      MAP_normal[49][i] = (0,1);
    }

  }
public:
	double 	stdDev(list<Particle> ParticleSet, int numofParticles, int stateSize);
	vector<double>	RandRRT(double*  map, int x_size, int y_size);
	void	nearestNeighborRRT(vector<double> qrand, vertex **nearestVertex);
	int 	extendRRT(vector<double> qrand, double*  map, int x_size, int y_size);
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
//Helper functions
double 	smallerAngle(double differ);
float 	distNorm(double* q, double* qrand1, int numofDOFs);
void 	print_vector(vector<double> vecToPrint, int numofDOFs);


#endif // TREE__
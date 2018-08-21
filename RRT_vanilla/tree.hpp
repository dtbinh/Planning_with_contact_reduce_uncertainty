/*
 * tree.hpp
 * Contains a class that defines tree structure
 */

#ifndef TREE__
#define TREE__

#include <list>
#include <vector>

#define PI 3.141592654
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;
typedef struct vertex{
    // double theta[5];
    vector<double> theta;
    // int vertexID; // Is this really needed?
    // double cost; //What is this cost?
    std::list<vertex*> children;
    // int parentID;
    vertex* parent;
} vertex;

typedef struct edge {
    double startID;
    double endID;
} edge;

class Tree
{

public:
	std::list<vertex*> vertices;
	double step_size;
	int numofDOFs;
	vertex* start = new vertex;
	vertex* goal = new vertex;
	float goalProbability;
	int maxconfigs;

	Tree(double* armstart_anglesV_rad,
	    double* armgoal_anglesV_rad,
	    int nDOFs){

		//Initializing start and goal
		start->parent = NULL;
		// start->vertexID = 1;
		goal->parent = NULL;
		step_size = 0.3;// PI*5/180; //% degrees increment
    	numofDOFs = nDOFs;
    	goalProbability = 0.90;
    	maxconfigs = 30000;
		//Blocks
		for(int i = 0 ; i < numofDOFs; i++){ 
			start->theta.push_back(*(armstart_anglesV_rad+i));
			goal->theta.push_back(*(armgoal_anglesV_rad+i));
    	}

    	vertices.push_back(start);
	}
public:
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
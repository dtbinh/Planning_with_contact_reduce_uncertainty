#include <math.h>
#include <ctime>
#include <chrono>
#include <random>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;

int main() 
{
	int stateSize = 3;
	const int numofParticles = 4;
    MatrixXd particleMatrix(stateSize,numofParticles);
    MatrixXd projectedTargetMatrix(stateSize,numofParticles);
    MatrixXd deltaMatrix(stateSize,numofParticles);
    VectorXd meanParticle(stateSize);
    MatrixXd covMat(stateSize,stateSize);
    vector<double> res;
    VectorXd normal(3);
    normal << 2, 2, 2;
    double CovCost;
    for (int i = 0; i < stateSize; i++){
	    for (int j = 0; j < numofParticles; j++)
	    {
	      particleMatrix(i,j) = i+j;
	    }
	}

	// VectorXi columns = 1;;
	// VectorXd extracted_cols = particleMatrix(Eigen::placeholders::all, columns);
	// cout << "\nExtracted cols\n" << extracted_cols << endl;
	// return 0;

	cout << "Particles\n" << particleMatrix << endl;
	meanParticle = particleMatrix.rowwise().mean();
	cout << "Mean particle\n" << meanParticle << endl;
	deltaMatrix = particleMatrix.colwise() - meanParticle;
	cout << "deltaMatrix\n" << deltaMatrix << endl;
	// stdDev = ((particleMatrix.colwise() - particleMatrix.rowwise().mean()).square().colwise().sum()/(numofParticles-1)).sqrt();
	covMat = (deltaMatrix * deltaMatrix.transpose())/numofParticles;
	cout << "covMat\n" << covMat << endl;
	CovCost = sqrt(covMat.trace());
	cout << "CovCost\n" << CovCost << endl;
	projectedTargetMatrix = (particleMatrix.cwiseProduct(normal.replicate<1,(numofParticles)>())).colwise().sum() ;
  	// projectedTargetMatrix = targetMatrix.colwise() - ( (targetMatrix.colwise() - nearestVertex->particleMatrix.colwise()).transpose()*normal)*normal;
  	// projectedTargetMatrix = ( (particleMatrix - deltaMatrix).transpose().rowwise()*normal)*normal;
  	cout << "normal.replicate<1,numofParticles>())\n" << normal.replicate<1,4>()	 << endl;
  	cout << "projectedTargetMatrix\n" << projectedTargetMatrix << endl;

  	vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd> > pathVec;
  	pathVec.push_back(particleMatrix);

  	VectorXd a(2);
  	VectorXd b(2);
  	VectorXd q(2);
  	a << 1,1;
  	b << 5,5;
  	double step = 2*1.414;
  	a = a + step*(b - a)/(b - a).norm();
  	cout << "q stepped\n" << a << endl;

	// cout << "Rand\n" << endl;
	// vector<double> rand = {7, 7, 7};
 //    for (int i = 0; i < rand.size(); i++) {
 //        cout << rand[i] << " ";
 //        res.push_back(rand[i] - meanParticle(i));
 //    }
 //    cout << "Res\n" << endl;
 //    for (int i = 0; i < rand.size(); i++) {
 //        cout << res[i] << " ";
        
 //    }

	
    return 0;
}
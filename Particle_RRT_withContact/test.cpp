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
	int numofParticles = 4;
    MatrixXd particleMatrix(stateSize,numofParticles);
    MatrixXd deltaMatrix(stateSize,numofParticles);
    VectorXd meanParticle(stateSize);
    MatrixXd covMat(stateSize,stateSize);
    vector<double> res;
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


	cout << "Rand\n" << endl;
	vector<double> rand = {7, 7, 7};
    for (int i = 0; i < rand.size(); i++) {
        cout << rand[i] << " ";
        res.push_back(rand[i] - meanParticle(i));
    }
    cout << "Res\n" << endl;
    for (int i = 0; i < rand.size(); i++) {
        cout << res[i] << " ";
        
    }

	
    return 0;
}
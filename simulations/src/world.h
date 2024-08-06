#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include input file and option
#include "setInput.h"

// include elastic Plate class
#include "elasticPlate.h"

// include time stepper
#include "timeStepper.h"

// include force
#include "inertialForce.h"
#include "externalGravityForce.h"
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "dampingForce.h"
#include "externalContactForce.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	void setPlateStepper();

	void updateTimeStep();

	int simulationRunning();

	int numStretchingPair();
	Vector2d getScaledCoordinate(int i, int j);
		
private:

	// physical parameters
	bool render;
	bool saveData;
	double deltaTime;
	double totalTime;
	double YoungM;
	double density;
	double Possion;
	double stol;
	double forceTol;
	double scaleRendering;
	int maxIter;
	Vector3d gVector;
	double viscosity;
	double rodRadius;
	double deltaLength;
	double epsilonV;
	double mu;
	double angleRight;
	double totalMass;

	double l1;
	double l2;
	double compressRatio;
	double h1;

	double stiffness;
    double dBar; 

	int Nstep;
	int timeStep;

	double characteristicForce;

	double currentTime;

	void plateBoundaryCondition();

	// Plate
	elasticPlate *plate;

	// stepper
	timeStepper *stepper;

	// force
	inertialForce *m_inertialForce;
	externalGravityForce *m_gravityForce;
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	dampingForce *m_dampingForce;
	externalContactForce *m_externalContactForce;

	void updateEachStep();

	std::vector<double> v_center;

	double totalEnergy;

	int ifHigher;
	int ifLower;

	double computeHeight(double x, double y);
	double computeHeight2(double x, double y);
};

#endif

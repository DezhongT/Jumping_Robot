#ifndef EXTERNALCONATCTFORCE_H
#define EXTERNALCONATCTFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar, double m_mu, double m_epsilonV);
	~externalContactForce();

	void computeFc();
	void computeJc();

    VectorXd ForceVec;

    int ifContact;

    double stiffness;
    double dBar;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    int ind;

    double dEnergydD;
    double d2EnergydD2;

    Vector2d f;

    double mu;
    double epsilonV;

    double dt;

    double fVelocity;
    Vector2d tK;

    Vector2d friction;
    Matrix2d frictionJacobian;
    Vector2d dfVelocity;
    Matrix2d dtK;

    Matrix2d Id3;
    Matrix2d IdG;
};

#endif

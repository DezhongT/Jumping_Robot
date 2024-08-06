#ifndef EXTERNALMAGNETICFORCE_H
#define EXTERNALMAGNETICFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalMagneticForce
{
public:
	externalMagneticForce(elasticPlate &m_plate, timeStepper &m_stepper, 
    Vector3d m_bAVector, Vector3d m_bRVector);
	~externalMagneticForce();

	void computeFm();
	void computeJm();

    VectorXd ForceVec;
	
private:
	elasticPlate *plate;
	timeStepper *stepper;

    Vector2d baVector_ref;
    Vector2d brVector_ref;

    Vector2d baVector;
    Vector2d brVector;

    double muZero;

    Vector2d m_current;
    Vector2d t_current;

    Vector2d m_start;
    Vector2d t_start;

    Matrix2d Id2;

    Matrix2d dtde;
    Matrix2d dmde;
    Vector2d dEde;

    Vector2d x1;
    Vector2d x2;

    Vector2d x1_start;
    Vector2d x2_start;

    double edge;

    VectorXd force;
    VectorXi localDOF;

    Vector2d M;

    int nv1;
    int nv2;
};

#endif

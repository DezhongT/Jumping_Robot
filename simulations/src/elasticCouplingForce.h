#ifndef ELASTICCOUPLINGFORCE_H
#define ELASTICCOUPLINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticCouplingForce
{
public:
	elasticCouplingForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticCouplingForce();

	void computeFcouple();
	void computeJcouple();

    void setFirstJacobian();

private:

	elasticPlate *plate;
    timeStepper *stepper;

    MatrixXd JacobLocal;
    VectorXd forceLocal;
    VectorXi localDOF;

    int ind1, ind2, ind3, ind4;
    Vector2d x1, x2, x3, x4;

    double l_k, phi0, k_local;

    VectorXd ListVec(double a1, double a2, double a3, 
    double a4, double a5, double a6, double a7, double a8);
    MatrixXd ListMat(VectorXd a1, VectorXd a2, VectorXd a3, 
    VectorXd a4, VectorXd a5, VectorXd a6, VectorXd a7, VectorXd a8);

    VectorXd computeCoupleForce(double x1, double y1, double x2, double y2, double x3, double y3, 
    double x4, double y4, double k, double phi0);
    MatrixXd computeCoupleJacobian(double x1, double y1, double x2, double y2, double x3, double y3, 
    double x4, double y4, double k, double phi0);

};

#endif

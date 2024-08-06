#include "elasticCouplingForce.h"

elasticCouplingForce::elasticCouplingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
    plate = &m_plate;
    stepper = &m_stepper;

    JacobLocal = MatrixXd::Zero(8,8);
    forceLocal = VectorXd::Zero(8);

    localDOF = VectorXi::Zero(8);
}

elasticCouplingForce::~elasticCouplingForce()
{
	;
}

void elasticCouplingForce::computeFcouple()
{
    for (int i = 0; i < plate->v_couplingPair.size(); i++)
    {
        ind1 = plate->v_couplingPair[i].nv_1;
        ind2 = plate->v_couplingPair[i].nv_2;
        ind3 = plate->v_couplingPair[i].nv_3;
        ind4 = plate->v_couplingPair[i].nv_4;

        x1 = plate->getVertex(ind1);
        x2 = plate->getVertex(ind2);
        x3 = plate->getVertex(ind3);
        x4 = plate->getVertex(ind4);

        l_k = plate->v_couplingPair[i].voroniLength;
        phi0 = plate->v_couplingPair[i].phi0;

        k_local = plate->v_couplingPair[i].k_local;

        forceLocal = - computeCoupleForce(x1(0), x1(1), x2(0), x2(1), x3(0), x3(1), x4(0), x4(1), k_local, phi0) / l_k;

        //cout << forceLocal.transpose() << endl;

        localDOF(0) = 2 * ind1 + 0;
        localDOF(1) = 2 * ind1 + 1;
        localDOF(2) = 2 * ind2 + 0;
        localDOF(3) = 2 * ind2 + 1;
        localDOF(4) = 2 * ind3 + 0;
        localDOF(5) = 2 * ind3 + 1;
        localDOF(6) = 2 * ind4 + 0;
        localDOF(7) = 2 * ind4 + 1;

        for (int k = 0; k < 6; k++)
        {
            stepper->addForce(localDOF(k), - forceLocal(k));
        }
    }
}

void elasticCouplingForce::computeJcouple()
{
    for (int i = 0; i < plate->v_couplingPair.size(); i++)
    {
        ind1 = plate->v_couplingPair[i].nv_1;
        ind2 = plate->v_couplingPair[i].nv_2;
        ind3 = plate->v_couplingPair[i].nv_3;
        ind4 = plate->v_couplingPair[i].nv_4;

        x1 = plate->getVertex(ind1);
        x2 = plate->getVertex(ind2);
        x3 = plate->getVertex(ind3);
        x4 = plate->getVertex(ind4);

        l_k = plate->v_couplingPair[i].voroniLength;
        phi0 = plate->v_couplingPair[i].phi0;

        k_local = plate->v_couplingPair[i].k_local;

        JacobLocal = computeCoupleJacobian(x1(0), x1(1), x2(0), x2(1), x3(0), x3(1), x4(0), x4(1), k_local, phi0) / l_k;

        localDOF(0) = 2 * ind1 + 0;
        localDOF(1) = 2 * ind1 + 1;
        localDOF(2) = 2 * ind2 + 0;
        localDOF(3) = 2 * ind2 + 1;
        localDOF(4) = 2 * ind3 + 0;
        localDOF(5) = 2 * ind3 + 1;
        localDOF(6) = 2 * ind4 + 0;
        localDOF(7) = 2 * ind4 + 1;

        for (int ii = 0; ii < 8; ii++)
        {
            for (int jj = 0; jj < 8; jj++)
            {
                stepper->addJacobian(localDOF(ii), localDOF(jj), JacobLocal(ii, jj));
            }
        }
    }
}

void elasticCouplingForce::setFirstJacobian()
{
    for (int i = 0; i < plate->v_couplingPair.size(); i++)
    {
        ind1 = plate->v_couplingPair[i].nv_1;
        ind2 = plate->v_couplingPair[i].nv_2;
        ind3 = plate->v_couplingPair[i].nv_3;
        ind4 = plate->v_couplingPair[i].nv_4;

        localDOF(0) = 2 * ind1 + 0;
        localDOF(1) = 2 * ind1 + 1;
        localDOF(2) = 2 * ind2 + 0;
        localDOF(3) = 2 * ind2 + 1;
        localDOF(4) = 2 * ind3 + 0;
        localDOF(5) = 2 * ind3 + 1;
        localDOF(6) = 2 * ind4 + 0;
        localDOF(7) = 2 * ind4 + 1;

        for (int ii = 0; ii < 8; ii++)
        {
            for (int jj = 0; jj < 8; jj++)
            {
                stepper->addJacobian(localDOF(ii), localDOF(jj), 1);
            }
        }
    }
}


VectorXd elasticCouplingForce::computeCoupleForce(double x1, double y1, double x2, double y2, double x3, double y3, 
    double x4, double y4, double k, double phi0)
{
    VectorXd vecResult;

    vecResult = ListVec((-1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
         (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
        (-x3 + x4)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
         (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
   (-1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
         (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        (pow(-y1 + y2,2)*(-y3 + y4))/
         (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
        (-y3 + y4)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
   (-1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
        (-x3 + x4)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
        ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
         (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
   (-1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
        (pow(-y1 + y2,2)*(-y3 + y4))/
         (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        (-y3 + y4)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
   (-1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
        ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
        (-x1 + x2)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
   (-1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
        ((-y1 + y2)*pow(-y3 + y4,2))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
        (-y1 + y2)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
   (-1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
        ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
        (-x1 + x2)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
   (-1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
        ((-y1 + y2)*pow(-y3 + y4,2))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
        (-y1 + y2)/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
      (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
    sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
        ((-y1 + y2)*(-y3 + y4))/
         (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
           sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)));

    return vecResult;
}

MatrixXd elasticCouplingForce::computeCoupleJacobian(double x1, double y1, double x2, double y2, double x3, double y3, 
    double x4, double y4, double k, double phi0)
{
    MatrixXd matResult;

    matResult = ListMat(ListVec((1.*k*pow((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*pow(-x1 + x2,3)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*pow(-x1 + x2,2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*pow(-x1 + x2,3)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*pow(-x1 + x2,2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))),
   ListVec((1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*pow(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*(-x3 + x4)*pow(-y1 + y2,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*pow(-y1 + y2,3)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*(-x3 + x4)*pow(-y1 + y2,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*pow(-y1 + y2,3)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))),
   ListVec((1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*pow(-x1 + x2,3)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*pow(-x1 + x2,2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*pow(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*pow(-x1 + x2,3)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*pow(-x1 + x2,2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))),
   ListVec((1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*(-x3 + x4)*pow(-y1 + y2,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*pow(-y1 + y2,3)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*pow(-x1 + x2,2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*(-x1 + x2)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*pow(-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow(-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*(-x3 + x4)*pow(-y1 + y2,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (3*pow(-y1 + y2,3)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),2.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (3*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))),
   ListVec((1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*pow(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*pow(-x3 + x4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*pow(-x3 + x4,2)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*pow(-x3 + x4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*pow(-x3 + x4,2)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))),
   ListVec((1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*pow(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*(-x3 + x4)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-y1 + y2)*pow(-y3 + y4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          ((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (3*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*(-x3 + x4)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-y1 + y2)*pow(-y3 + y4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          ((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (3*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))),
   ListVec((1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*pow(-x3 + x4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x3 + x4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-x1 + x2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*pow(-x3 + x4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*pow(-x3 + x4,2)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*pow(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*pow(-x3 + x4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*pow(-x3 + x4,2)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))),
   ListVec((1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        ((pow(-x1 + x2,2)*(-x3 + x4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) + 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-((pow(-x1 + x2,2)*(-x3 + x4))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) + 
          (-x3 + x4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) - 
          ((-x1 + x2)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((pow(-x1 + x2,2)*(-x3 + x4)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y1 + y2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y1 + y2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y1 + y2))/
             (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))) - 
          (pow(-y1 + y2,2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          (-y3 + y4)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (pow(-y1 + y2,2)*pow(-y3 + y4,2))/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y3 + y4,2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          pow(-y1 + y2,2)/
           (pow(pow(-x1 + x2,2) + pow(-y1 + y2,2),1.5)*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          1/(sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*pow(-x3 + x4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((-3*(-x1 + x2)*(-x3 + x4)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          (3*(-y1 + y2)*pow(-y3 + y4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          ((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (3*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2)))))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*(-(((-x1 + x2)*pow(-x3 + x4,2))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-x3 + x4)*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-x1 + x2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*pow(-x3 + x4,2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-x3 + x4)*(-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          ((-x3 + x4)*(-y1 + y2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          ((-x1 + x2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)),
    (1.*k*pow(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))/
      (1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)) - 
     (1.*k*pow(-(((-x1 + x2)*(-x3 + x4)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5))) - 
          ((-y1 + y2)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) + 
          (-y1 + y2)/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2)*
        (((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      pow(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2),1.5) - 
     (1.*k*((3*(-x1 + x2)*(-x3 + x4)*pow(-y3 + y4,2))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) + 
          (3*(-y1 + y2)*pow(-y3 + y4,3))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),2.5)) - 
          ((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)) - 
          (3*(-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             pow(pow(-x3 + x4,2) + pow(-y3 + y4,2),1.5)))*
        (-phi0 + acos(((-x1 + x2)*(-x3 + x4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
            ((-y1 + y2)*(-y3 + y4))/
             (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
               sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))))))/
      sqrt(1 - pow(((-x1 + x2)*(-x3 + x4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))) + 
          ((-y1 + y2)*(-y3 + y4))/
           (sqrt(pow(-x1 + x2,2) + pow(-y1 + y2,2))*
             sqrt(pow(-x3 + x4,2) + pow(-y3 + y4,2))),2))));

    return matResult;
}

VectorXd elasticCouplingForce::ListVec(double a1, double a2, double a3, 
    double a4, double a5, double a6, double a7, double a8)
{
    VectorXd vecResult;

    vecResult.setZero(8, 1);

    vecResult(0) = a1;
    vecResult(1) = a2;
    vecResult(2) = a3;
    vecResult(3) = a4;
    vecResult(4) = a5;
    vecResult(5) = a6;
    vecResult(6) = a7;
    vecResult(7) = a8;

    return vecResult;
}

MatrixXd elasticCouplingForce::ListMat(VectorXd a1, VectorXd a2, VectorXd a3, 
    VectorXd a4, VectorXd a5, VectorXd a6, VectorXd a7, VectorXd a8)
{
    MatrixXd matResult;

    matResult.setZero(8, 8);

    matResult.col(0) = a1;
    matResult.col(1) = a2;
    matResult.col(2) = a3;
    matResult.col(3) = a4;
    matResult.col(4) = a5;
    matResult.col(5) = a6;
    matResult.col(6) = a7;
    matResult.col(7) = a8;

    return matResult;
}
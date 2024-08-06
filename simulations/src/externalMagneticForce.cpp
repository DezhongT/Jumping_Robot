#include "externalMagneticForce.h"

externalMagneticForce::externalMagneticForce(elasticPlate &m_plate, timeStepper &m_stepper, 
	Vector3d m_bAVector, Vector3d m_bRVector)
{
	plate = &m_plate;
	stepper = &m_stepper;

	baVector_ref(0) = m_bAVector(0);
	baVector_ref(1) = m_bAVector(1);

	brVector_ref(0) = m_bRVector(0);
	brVector_ref(1) = m_bRVector(1);

	muZero = 1.0;

	Id2<<1,0,
         0,1;

    force.setZero(4, 1);
    localDOF.setZero(4, 1);

    ForceVec = VectorXd::Zero(plate->ndof);
}

externalMagneticForce::~externalMagneticForce()
{
	;
}

void externalMagneticForce::computeFm()
{
	/*
	ForceVec = VectorXd::Zero(plate->ndof);

	for (int i = 0; i < plate->edgeNum; i++)
	{
		nv1 = plate->v_edgeElement[i].nv_1;
		nv2 = plate->v_edgeElement[i].nv_2;

		x1 = plate->getVertexOld(nv1); 
		x2 = plate->getVertexOld(nv2); 

		brVector.setZero(2, 1);

		if (i < plate->curveNv / 3)
		{
			brVector = brVector_ref;
		}

		if (i > plate->curveNv / 3 && i < plate->curveNv)
		{
			brVector = - brVector_ref;
		}
	
		if (i > plate->frameNvS && i < plate->frameNvS + (plate->frameNvE - plate->frameNvS) / 2)
		{
			brVector = - brVector / 10;
		}

		if (i > plate->frameNvE && i > plate->frameNvS + (plate->frameNvE - plate->frameNvS) / 2)
		{
			brVector = brVector / 10;
		}

		baVector(0) = baVector_ref(0);
		baVector(1) = baVector_ref(1);

		t_current = (x2 - x1) / (x2 - x1).norm();
		m_current(0) = - t_current(1);
		m_current(1) = t_current(0);

		edge = (x2 - x1).norm();

		x1_start = plate->getVertexStart(nv1); 
		x2_start = plate->getVertexStart(nv2); 

		t_start = (x2_start - x1_start) / (x2_start - x1_start).norm();
		m_start(0) = - t_start(1);
		m_start(1) = t_start(0);

		dtde = (Id2 - t_current * t_current.transpose()) / edge;
		dmde = - (t_current * m_current.transpose()) / edge;

		M = m_start.dot(brVector) * m_current + t_start.dot(brVector) * t_current;

		dEde = (m_start.dot(brVector) * dmde + t_start.dot(brVector) * dtde).transpose() * baVector;

		force.segment(0, 2) = - dEde;
		force.segment(2, 2) = dEde;

		force = - force * edge * plate->crossSectionalArea / muZero;

		localDOF(0) = 2 * nv1 + 0;
		localDOF(1) = 2 * nv1 + 1;
		localDOF(2) = 2 * nv2 + 0;
		localDOF(3) = 2 * nv2 + 1;

		for (int i = 0; i < 4; i++)
		{
			stepper->addForce(localDOF(i), - force(i));
		}

	}
	*/

	//int index1 = plate->curveNv / 2;
	//cout << index1 << " " << - brVector_ref(0) << endl;
	//stepper->addForce(2 * index1 + 1, - brVector_ref(0));

	//int index2 = (plate->frameNvE + plate->frameNvS) / 2;
	//stepper->addForce(2 * index2 + 1, - baVector_ref(0));
}

void externalMagneticForce::computeJm()
{
	;
}

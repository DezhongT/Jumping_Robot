#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar, double m_mu, double m_epsilonV)
{
	plate = &m_plate;
	stepper = &m_stepper;

    stiffness = m_stiffness;
    dBar = m_dBar;

    mu = m_mu;

    epsilonV = m_epsilonV;

    dt = plate->dt;

    ForceVec = VectorXd::Zero(plate->ndof);

    Id3<<1,0,
         0,1;

    IdG<<1,0,
    	 0,0;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{
	ifContact = 0;

	ForceVec = VectorXd::Zero(plate->ndof);

	for(int i = 0; i < plate->nv; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);

		double d = xCurrent(1);

		if (d <= dBar)
		{
			//cout << "height :" << d << endl;

			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			stepper->addForce(2 * i + 1, stiffness * dEnergydD);

			ifContact = 2;
		}

		
		xCurrent = plate->getVertexOld(i);

		d = xCurrent(1);

		if (d <= dBar)
		{
			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			Vector2d uCurrent = plate->getVelocityOld(i);

			Vector2d vTangent;

			vTangent(0) = uCurrent(0);
			vTangent(1) = 0.0;
			//vTangent(2) = 0.0;

			if ( vTangent.norm() >= epsilonV )
			{
				fVelocity = 1.0;
				tK = vTangent / vTangent.norm();
			}
			else
			{
				fVelocity = - ( vTangent.norm() * vTangent.norm() ) / (epsilonV * epsilonV) + 2 * vTangent.norm() / epsilonV;
				tK = vTangent / (vTangent.norm() + 1e-15);
			}

			friction = - mu * abs(stiffness * dEnergydD) * fVelocity * tK;

			for (int kk = 0; kk < 2; kk++)
			{
				stepper->addForce(2 * i + kk, - friction(kk));
			}
		}

	}

	//stepper->subtractForceVector(-ForceVec);
}

void externalContactForce::computeJc()
{
	for(int i = 0; i < plate->nv; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);

		double d = xCurrent(1);

		if (d <= dBar)
		{
			d2EnergydD2 = - 2 * log(d / dBar) - 2 * (d - dBar) / d - 2 * (d - dBar) / d + (d - dBar) * (d - dBar) / (d * d);

			stepper->addJacobian(2 * i + 1, 2 * i + 1, stiffness * d2EnergydD2);
		}

		//xCurrent = plate->getVertexOld(i);

		//d = xCurrent(2);

		/*

		if (d <= dBar)
		{
			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			Vector3d uCurrent = plate->getVelocity(i);

			Vector3d vTangent;

			vTangent(0) = uCurrent(0);
			vTangent(1) = uCurrent(1);
			vTangent(2) = 0.0;

			if ( vTangent.norm() >= epsilonV )
			{
				fVelocity = 1.0;
				tK = vTangent / vTangent.norm();
			}
			else
			{
				fVelocity = - ( vTangent.norm() * vTangent.norm() ) / (epsilonV * epsilonV) + 2 * vTangent.norm() / epsilonV;
				tK = vTangent / (vTangent.norm() + 1e-15);
			}

			if ( vTangent.norm() >= epsilonV )
			{
				dfVelocity(0) = 0.0;
				dfVelocity(1) = 0.0;
				dfVelocity(2) = 0.0;

				dtK = (Id3 - tK * tK.transpose() ) / vTangent.norm();
			}
			else
			{
				dfVelocity = ( - 2 * vTangent.norm() / (epsilonV * epsilonV) + 2 / epsilonV ) * tK;

				dtK = (Id3 - tK * tK.transpose() ) / (vTangent.norm() + 1e-15);
			}

			if (dEnergydD > 0.0)
			{
				d2EnergydD2 = d2EnergydD2;
			}
			else
			{
				d2EnergydD2 = - d2EnergydD2;
			}

			//friction = mu * abs(stiffness * dEnergydD) * fVelocity * tK;

			Vector3d d2EnergydD2Local;

			d2EnergydD2Local(0) = 0.0;
			d2EnergydD2Local(1) = 0.0;
			d2EnergydD2Local(2) = d2EnergydD2;

			frictionJacobian = mu * stiffness * fVelocity * d2EnergydD2Local * tK.transpose() + mu * abs(stiffness * dEnergydD) * ( dfVelocity * tK.transpose() + fVelocity * dtK) * IdG / plate->dt;

			for (int jj = 0; jj < 3; jj++)
			{
				for (int kk = 0; kk < 3; kk++)
				{
					stepper->addJacobian(4 * i + jj, 4 * i + kk, frictionJacobian(jj, kk) );
				}
			}
		}
		*/
	}
}

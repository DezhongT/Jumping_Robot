#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");			
	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");
	rodRadius = m_inputData.GetScalarOpt("rodRadius");
	Possion = m_inputData.GetScalarOpt("Possion");
	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");
	deltaLength = m_inputData.GetScalarOpt("deltaLength");
	stiffness = m_inputData.GetScalarOpt("stiffness");
	dBar = m_inputData.GetScalarOpt("dBar");
	totalMass = m_inputData.GetScalarOpt("totalMass");

	l1 = m_inputData.GetScalarOpt("l1");
	l2 = m_inputData.GetScalarOpt("l2");
	compressRatio = m_inputData.GetScalarOpt("compressRatio");
	h1 = m_inputData.GetScalarOpt("h1");
	angleRight = m_inputData.GetScalarOpt("angleRight");

	//cout << computeHeight(1 - compressRatio, angleRight) << endl;
	//cout << computeHeight2(1 - compressRatio, angleRight) << endl;

	h1 = - 1.01 * l1 * ( computeHeight2(1 - compressRatio, angleRight) + computeHeight2(1 - compressRatio, angleRight) ) / 2 ;

	//cout << compressRatio << " " << h1 << endl;

	mu = m_inputData.GetScalarOpt("mu");
	epsilonV = m_inputData.GetScalarOpt("epsilonV");


	v_center.clear();

	totalEnergy = 0.0;

	ifHigher = 1;
	ifLower = 1;
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDiscreteNet";
    name << "_l2_" << l2;
    name << "_compressRatio_" << compressRatio;
    name << "_angleRight_" << angleRight;
    name << "_mu_" << mu;
    name << "_mass_" << totalMass;
    //name << "_h1_" << h1;
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if ( timeStep % 10 != 0)
	{
		return;
	}

	if (timeStep == Nstep)
	{
		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);

			//outfile << xCurrent(0) << " " << xCurrent(1) << endl;
		}

		for (int i = 0; i < plate->edge.size(); i++)
		{
			Vector2i edgeCurrent = plate->edge[i];

			//cout << edgeCurrent(0) << " " << edgeCurrent(1) << endl;
		}
	}

	if (currentTime > 5.0)
	{
		if (m_externalContactForce->ifContact > 1)
		{
			timeStep = Nstep;
		}
	}

	if (currentTime > 3.9)
	{
		Vector2d xCenter;

		xCenter(0) = 0.0;
		xCenter(1) = 0.0;

		//int indexOut = 0;

		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);

			double mass = plate->massArray(2 * i);

			//outfile << currentTime - 4.0 << " " << xCurrent(0) << " " << xCurrent(1) << endl;

			xCenter = xCenter + xCurrent * mass;
		}

		xCenter = xCenter / (plate->totalM);

		outfile << currentTime - 4.0 << " " << xCenter(0) << " " << xCenter(1) << endl;

		/*

		m_bendingForce->computeMoment();

		if (m_bendingForce->boundaryMoment > 0.0)
		{
			totalEnergy = totalEnergy + deltaTime * 0.1 * m_bendingForce->boundaryMoment;
		}

		outfile << totalEnergy << " " << xCenter(1) << endl;

		*/

		//m_bendingForce->computeMoment();

		//outfile << currentTime << " " << atan(plate->v_bendingElement[plate->bendingNum - 1].nBar / 2) * 2 << " " << m_bendingForce->boundaryMoment << endl;



		/*

		v_center.push_back(xCenter(1));
		if (v_center.size() > 10 && currentTime > 6.2)
		{
			if ( v_center[v_center.size() - 1] < v_center[v_center.size() - 5])
			{
				if (m_externalContactForce->ifContact < 1 && abs(plate->v_bendingElement[plate->bendingNum - 1].nBar) > 0.5)
				{
					timeStep = Nstep;
				}
			}
		}

		*/
	}

	if (timeStep == Nstep)
	{
		//outfile << v_center[v_center.size() - 1] << endl;
	}

}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, rodRadius, Possion, deltaTime, deltaLength, l1, l2, compressRatio, h1, totalMass);

	plateBoundaryCondition();

	plate->setup();

	stepper = new timeStepper(*plate);

	// set up force
	m_inertialForce = new inertialForce(*plate, *stepper);
	m_gravityForce = new externalGravityForce(*plate, *stepper, gVector);
	m_stretchForce = new elasticStretchingForce(*plate, *stepper);
	m_bendingForce = new elasticBendingForce(*plate, *stepper);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosity);
	m_externalContactForce = new externalContactForce(*plate, *stepper, stiffness, dBar, mu, epsilonV);
	
	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_stretchForce->setFirstJacobian();
	m_bendingForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();
	
	stepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;
}

void world::plateBoundaryCondition()
{
	plate->setVertexBoundaryCondition(plate->getVertex(0), 0);
	plate->setVertexBoundaryCondition(plate->getVertex(1), 1);

	plate->setVertexBoundaryCondition(plate->getVertex(plate->curveNv-2), plate->curveNv-2);
	plate->setVertexBoundaryCondition(plate->getVertex(plate->curveNv-1), plate->curveNv-1);

	for (int i = plate->curveNv; i < plate->nv; i++)
	{
		plate->setVertexBoundaryCondition(plate->getVertex(i), i);
	}
}

void world::updateTimeStep()
{
	bool goodSolved = false;

	while (goodSolved == false)
	{
		// Start with a trial solution for our solution x
		plate->updateGuess(); // x = x0 + u * dt

		updateEachStep();

		goodSolved = true;
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " ";
	}

	currentTime += deltaTime;
		
	timeStep++;
}

void world::updateEachStep()
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	double beamLength = ( plate->getVertex(plate->curveNv-1) - plate->getVertex(0) ).norm();

	if (currentTime > 0.05 && currentTime < 1.0 && beamLength > plate->l1 - 1.99 * plate->deltaL)
	{
		Vector2d x1 = plate->getVertex(0);
		x1(0) = x1(0) + deltaTime * 0.01;
		plate->setVertexBoundaryCondition(x1, 0);

		Vector2d x2 = plate->getVertex(1);
		x2(0) = x2(0) + deltaTime * 0.01;
		plate->setVertexBoundaryCondition(x2, 1);

		Vector2d x3 = plate->getVertex(plate->curveNv-2);
		x3(0) = x3(0) - deltaTime * 0.01;
		plate->setVertexBoundaryCondition(x3, plate->curveNv-2);

		Vector2d x4 = plate->getVertex(plate->curveNv-1);
		x4(0) = x4(0) - deltaTime * 0.01;
		plate->setVertexBoundaryCondition(x4, plate->curveNv-1);
	}

	if (timeStep == 1000)
	{
		plate->clearMap();

		plate->addElement();

		for (int i = plate->curveNv; i < plate->nv; i++)
		{
			plate->setVertexBoundaryCondition(plate->getVertex(i), i);
		}

		plate->reSetupMap();
		stepper->resetMap();

		// set up first jacobian
		m_inertialForce->setFirstJacobian();
		m_stretchForce->setFirstJacobian();
		m_bendingForce->setFirstJacobian();
		m_dampingForce->setFirstJacobian();

		stepper->first_time_PARDISO_setup();
	}


	if (timeStep == 2000)
	{
		plate->clearMap();

		plate->reSetupMap();
		stepper->resetMap();

		// set up first jacobian
		m_inertialForce->setFirstJacobian();
		m_stretchForce->setFirstJacobian();
		m_bendingForce->setFirstJacobian();
		m_dampingForce->setFirstJacobian();

		stepper->first_time_PARDISO_setup();

		m_gravityForce->gVector(0) = 0.0;
		m_gravityForce->gVector(1) = -10.0;
		m_gravityForce->setGravity();

		/*

		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);

			cout << xCurrent(0) << " " << xCurrent(1) << endl;
		}

		*/
	}


	if (timeStep == 3000)
	{
		Vector2d xMid = plate->getVertex( int(plate->curveNv - 1) / 2 );
		Vector2d x0 = plate->getVertex( 0 );

		if (xMid(1) < plate->l2)
		{
			ifHigher = 1;
		}
		else
		{
			ifHigher = 0;
		}

		if ( (xMid(1) - x0(1)) > x0(1) )
		{
			ifLower = 1;
		}
		else
		{
			ifLower = 0;
		}


		for (int i = 0; i < plate->curveNv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);

			if (xCurrent(1) < dBar)
			{
				ifLower = 0;
			}
		}

	}



	if (currentTime > 1.0 && currentTime < 3.0)
	{
		if (ifHigher == 0)
		{
			timeStep = Nstep;
		}

		if (ifLower == 0)
		{
			timeStep = Nstep;
		}

		if (angleRight > 0.0)
		{
			Vector2d xRight1 = plate->getVertex(plate->curveNv - 2);
			Vector2d xRight2 = plate->getVertex(plate->curveNv - 1);
			double angleRightCurrent = atan2(xRight2(1) - xRight1(1), xRight2(0) - xRight1(0));

			if ( abs(angleRightCurrent) < abs(angleRight) )
			{
				plate->v_bendingElement[plate->bendingNum - 1].nBar = plate->v_bendingElement[plate->bendingNum - 1].nBar + 2 * tan(deltaTime * 0.1 / 2);
			}
		}

		if (angleRight < 0.0)
		{
			Vector2d xRight1 = plate->getVertex(plate->curveNv - 2);
			Vector2d xRight2 = plate->getVertex(plate->curveNv - 1);
			double angleRightCurrent = atan2(xRight2(1) - xRight1(1), xRight2(0) - xRight1(0));

			if ( abs(angleRightCurrent) < abs(angleRight) )
			{
				plate->v_bendingElement[plate->bendingNum - 1].nBar = plate->v_bendingElement[plate->bendingNum - 1].nBar - 2 * tan(deltaTime * 0.1 / 2);
			}
		}
	}

	if (currentTime > 4.0 && currentTime < 4.075)
	{
		deltaTime = 1e-5;
		plate->dt = 1e-5;

		plate->v_bendingElement[plate->bendingNum - 1].nBar = plate->v_bendingElement[plate->bendingNum - 1].nBar - 2 * tan(deltaTime * 20.0 / 2);
		plate->v_bendingElement[plate->bendingNum - 3].nBar = plate->v_bendingElement[plate->bendingNum - 3].nBar + 2 * tan(deltaTime * 20.0 / 2);
	}


	if (currentTime > 4.0)
	{
		m_dampingForce->viscosity = 0.01;
	}

	if (currentTime > 4.1)
	{
		deltaTime = 1e-3;
		plate->dt = 1e-3;
	}




	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		m_gravityForce->computeFg();
		m_stretchForce->computeFs();
		m_bendingForce->computeFb();
		m_dampingForce->computeFd();
		m_externalContactForce->computeFc();
		
		normf = stepper->GlobalForceVec.norm();

		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			m_gravityForce->computeJg();
			m_stretchForce->computeJs();
			m_bendingForce->computeJb();
			m_dampingForce->computeJd();
			m_externalContactForce->computeJc();
			
			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			timeStep = Nstep;
			break;
		}
	}

	if (render)
	{
		cout << "iter " << iter << endl;
	}
}

int world::simulationRunning()
{
	if (currentTime < 10.0 && timeStep < Nstep) 
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

Vector2d world::getScaledCoordinate(int i, int j)
{
	Vector2d xCurrent;
	
	if (j == 0)
	{
		xCurrent = plate->v_edgeElement[i].x_1 * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = plate->v_edgeElement[i].x_2 * scaleRendering;
	}

	return xCurrent;
}

int world::numStretchingPair()
{
	return plate->edgeNum;
}



double world::computeHeight(double x, double y)
{
	/*

	double p00 =    -0.02461;
    double p10 =     -0.4186;
    double p01 =    -0.01057;
    double p20 =      0.9979;
    double p11 =      0.0638;
    double p02 =    -0.03302;
    double p30 =      -1.485;
    double p21 =      -0.203;
    double p12 =      0.1174;
    double p03 =    -0.03102;
    double p40 =       1.033;
    double p31 =      0.1756;
    double p22 =     -0.1428;
    double p13 =     0.07214;
    double p04 =   -0.002648;

	double height;

	height = p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p30*x*x*x + p21*x*x*y + p12*x*y*y + p03*y*y*y + p40*pow(x,4) + p31*pow(x,3)*y + p22*pow(x,2)*pow(y,2) + p13*x*pow(y,3) + p04*pow(y,4);

	height = - height;

	//cout << height << endl;

    return height;

    */

    double p1 = 1.18;
    double p2 = -1.065;
    double p3 = -0.1026;

    double height;

    height = p1*x*x + p2*x + p3;

    return height;
}


double world::computeHeight2(double x, double y)
{
	/*
	double p00 =    -0.07875;
    double p10 =      -1.485;
    double p01 =     0.04214;
    double p20 =       3.615 ;
    double p11 =     -0.1158 ;
    double p02 =     0.09813  ;
    double p30 =      -5.376  ;
    double p21 =      0.3602  ;
    double p12 =      -0.339  ;
    double p03 =      0.0437  ;
    double p40 =       3.413  ;
    double p31 =     -0.3174  ;
    double p22 =      0.3736  ;
    double p13 =    -0.09802  ;
    double p04 =      0.0141  ;

	double height;

	height = p00 + p10*x + p01*y + p20*x*x + p11*x*y + p02*y*y + p30*x*x*x + p21*x*x*y + p12*x*y*y + p03*y*y*y + p40*pow(x,4) + p31*pow(x,3)*y + p22*pow(x,2)*pow(y,2) + p13*x*pow(y,3) + p04*pow(y,4);

	height = - height;

	//cout << height << endl;

    return height;
    */

    double p1 = 0.3707;
    double p2 = - 0.3156;
    double p3 = - 0.02951;
    

    double height;

    height = p1*x*x + p2*x + p3;

    return height;

}
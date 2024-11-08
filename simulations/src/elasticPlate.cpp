#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_radius, 
		double m_Possion, double m_dt, double m_deltaLength,
		 double m_l1, double m_l2, double m_compressRatio, double m_h1, double m_mass)
{
	YoungM = m_YoungM;
	density = m_density;
	radius = m_radius;
	Possion = m_Possion;
	dt = m_dt;

	totalMass = m_mass / 5e-3;

	l1 = m_l1 * m_compressRatio;
    l2 = m_l2;
    compressRatio = m_compressRatio;
    beaml = m_l1;
    h1 = m_h1;

    deltaL = m_deltaLength;

	EA = YoungM * radius;
	EI = YoungM * radius * radius * radius / 12;
	crossSectionalArea = radius;

	setupGeometry();

	ndof = 2 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	for (int i = 0; i < nv; i++)
	{
		x(2 * i + 0) = v_nodes[i](0);
		x(2 * i + 1) = v_nodes[i](1) + 1e-3;
	}
	x0 = x;
	xSta = x;

	computeEdge();
	computeBending();

	setupMass();

	//set up constraint map
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{
	massArray = VectorXd::Zero(ndof);

	double deltaMass;

	int index1;
	int index2;

	double totalMass1 = 0.0;

	for (int i = 0; i < edgeNum; i++)
	{
		deltaMass = radius * density * v_edgeElement[i].refLength / 2;

		index1 = v_edgeElement[i].nv_1;
		index2 = v_edgeElement[i].nv_2;

		massArray(2 * index1 + 0) = massArray(2 * index1 + 0) + deltaMass;
		massArray(2 * index1 + 1) = massArray(2 * index1 + 1) + deltaMass;
	
		massArray(2 * index2 + 0) = massArray(2 * index2 + 0) + deltaMass;
		massArray(2 * index2 + 1) = massArray(2 * index2 + 1) + deltaMass;

		totalMass1 = totalMass1 + 2 * deltaMass;
	}

	deltaMass = (totalMass - totalMass1) / (nv - curveNv - 1);

	for (int i = curveNv; i < nv; i++)
	{
		massArray(2 * i + 0) = massArray(2 * i + 0) + deltaMass;
		massArray(2 * i + 1) = massArray(2 * i + 1) + deltaMass;
	}

	totalM = 0.0;

	for (int i = 0; i < nv; i++)
	{
		totalM = totalM + massArray(2 * i + 0) + massArray(2 * i + 1);
	}

	totalM = totalM / 2;

	//cout << massArray << endl;

	//cout << "Total mass : " << totalM * 5e-3 << endl;

}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector2d position, int k)
{
	isConstrained[2 * k + 0] = 1;
	isConstrained[2 * k + 1] = 1;
	
	// Store in the constrained dof vector
	x(2 * k + 0) = position(0);
	x(2 * k + 1) = position(1);
}

Vector2d elasticPlate::getVertex(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x(2 * i + 0);
	xCurrent(1) = x(2 * i + 1);
	
	return xCurrent;
}

Vector2d elasticPlate::getVertexStart(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = xSta(2 * i + 0);
	xCurrent(1) = xSta(2 * i + 1);
	
	return xCurrent;
}

Vector2d elasticPlate::getVertexOld(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x0(2 * i + 0);
	xCurrent(1) = x0(2 * i + 1);

	return xCurrent;
}

Vector2d elasticPlate::getVelocity(int i)
{
	Vector2d uCurrent;

	uCurrent(0) = ( x(2 * i + 0) - x0(2 * i + 0) ) / dt;
	uCurrent(1) = ( x(2 * i + 1) - x0(2 * i + 1) ) / dt;
	
	return uCurrent;
}

Vector2d elasticPlate::getVelocityOld(int i)
{
	Vector2d uCurrent;

	uCurrent(0) = u(2 * i + 0);
	uCurrent(1) = u(2 * i + 1);
	
	return uCurrent;
}

void elasticPlate::updateEdgePair()
{
	for (int i = 0; i < edgeNum; i++)
	{
		v_edgeElement[i].x_1 = getVertex(v_edgeElement[i].nv_1);
		v_edgeElement[i].x_2 = getVertex(v_edgeElement[i].nv_2);
		v_edgeElement[i].edgeLength = (v_edgeElement[i].x_1 - v_edgeElement[i].x_2).norm();
	}
}

void elasticPlate::updateBendingPair()
{
	for (int i = 0; i < bendingNum; i++)
	{
		v_bendingElement[i].x_1 = getVertex(v_bendingElement[i].nv_1);
		v_bendingElement[i].x_2 = getVertex(v_bendingElement[i].nv_2);
		v_bendingElement[i].x_3 = getVertex(v_bendingElement[i].nv_3);

		v_bendingElement[i].e_1 = v_bendingElement[i].x_2 - v_bendingElement[i].x_1;
		v_bendingElement[i].e_2 = v_bendingElement[i].x_3 - v_bendingElement[i].x_2;

		v_bendingElement[i].norm_1 =  v_bendingElement[i].e_1.norm();
		v_bendingElement[i].norm_2 =  v_bendingElement[i].e_2.norm();

		v_bendingElement[i].t_1 = v_bendingElement[i].e_1 / v_bendingElement[i].norm_1;
		v_bendingElement[i].t_2 = v_bendingElement[i].e_2 / v_bendingElement[i].norm_2;
	}
}

void elasticPlate::prepareForIteration()
{
	updateEdgePair();
	updateBendingPair();
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
	}
}

void elasticPlate::setupGeometry()
{
	v_nodes.clear();
    edge.clear();
    bending.clear();
    constraint.clear();

    int nv1 = beaml / deltaL;
    int nv2 = l1 / deltaL;
    int nv3 = l2 / deltaL;

    double deltaL1 = beaml / nv1;
    double deltaL2 = l1 / nv2;
    double deltaL3 = l2 / nv3;

    // setup for curve beam
    for (int i = 0; i < nv1 + 1; i++)
    {
    	Vector2d xCurrent;

    	xCurrent(0) = i * deltaL1 - beaml / 2;
    	xCurrent(1) = h1;

    	v_nodes.push_back(xCurrent);
    }
    curveNv = v_nodes.size();

    // setup for frame
    for (int i = 0; i < nv3; i++)
    {
    	Vector2d xCurrent;

    	xCurrent(0) = - l1 / 2;
    	xCurrent(1) = i * deltaL3;

    	v_nodes.push_back(xCurrent);
    }

    frameNvS = v_nodes.size();

    for (int i = 0; i < nv2; i++)
    {
    	Vector2d xCurrent;

    	xCurrent(0) = i * deltaL2 - l1 / 2;
    	xCurrent(1) = l2;

    	if (i == 0)
    	{
    		xCurrent(0) = xCurrent(0) + deltaL2 / 2;
    		xCurrent(1) = xCurrent(1) - deltaL3 / 2;
    	}

    	v_nodes.push_back(xCurrent);
    }

    frameNvE = v_nodes.size();

    for (int i = 0; i < nv3 + 1; i++)
    {
    	Vector2d xCurrent;

    	xCurrent(0) = l1 / 2;
    	xCurrent(1) = l2 - i * deltaL3;

    	if (i == 0)
    	{
    		xCurrent(0) = xCurrent(0) - deltaL2 / 2;
    		xCurrent(1) = xCurrent(1) - deltaL3 / 2;
    	}

    	v_nodes.push_back(xCurrent);
    }

    nv = v_nodes.size();

    /*


    cout << " node : " << endl;
    for (int i = 0; i < v_nodes.size(); i++)
    {
    	Vector2d xCurrent = v_nodes[i];

    	cout << i << " " << xCurrent(0) << " " << xCurrent(1) << endl;
    }

    cout << curveNv << " " << frameNvS << " " << frameNvE << endl;

    */


    // setup edge
    for (int i = 0; i < curveNv - 1; i++)
    {
    	Vector2i edgeCurrent;

    	edgeCurrent(0) = i;
    	edgeCurrent(1) = i + 1;

    	edge.push_back(edgeCurrent);
    }

    for (int i = curveNv; i < nv - 1; i++)
    {
    	Vector2i edgeCurrent;

    	edgeCurrent(0) = i;
    	edgeCurrent(1) = i + 1;

    	edge.push_back(edgeCurrent);
    }

    edgeNum = edge.size();

    /*

    cout << " edge : " << endl;
    for (int i = 0; i < edge.size(); i++)
    {
    	Vector2i edgeCurrent = edge[i];

    	cout << edgeCurrent(0) << " " << edgeCurrent(1) << endl;
    }

    */


    // setup bending
    for (int i = 0; i < curveNv - 2; i++)
    {
    	Vector3i bendingCurrent;

    	bendingCurrent(0) = i;
    	bendingCurrent(1) = i + 1;
    	bendingCurrent(2) = i + 2;

    	bending.push_back(bendingCurrent);
    }

    for (int i = curveNv; i < nv - 2; i++)
    {
    	Vector3i bendingCurrent;

    	bendingCurrent(0) = i;
    	bendingCurrent(1) = i + 1;
    	bendingCurrent(2) = i + 2;

    	bending.push_back(bendingCurrent);
    }

    bendingNum = bending.size();

    /*

    cout << " bending : " << endl;
    for (int i = 0; i < bendingNum; i++)
    {
    	Vector3i bendingCurrent = bending[i];

    	cout << bendingCurrent(0) << " " << bendingCurrent(1) << " " << bendingCurrent(2) << endl;
    }

    */
}

void elasticPlate::computeEdge()
{

	edgeNum = 0;
	v_edgeElement.clear();

	for (int i = 0; i < edge.size(); i++)
	{
		Vector2i edgeCurrent = edge[i];

		edgeElement m_edgeElement;

		m_edgeElement.nv_1 = edgeCurrent(0);
		m_edgeElement.nv_2 = edgeCurrent(1);

		m_edgeElement.x_1 = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2 = getVertex(m_edgeElement.nv_2);

		m_edgeElement.refLength = (m_edgeElement.x_2- m_edgeElement.x_1).norm();
		m_edgeElement.edgeLength = m_edgeElement.refLength;

		m_edgeElement.arrayNum = VectorXi::Zero(4);

		m_edgeElement.arrayNum(0) = 2 * m_edgeElement.nv_1 + 0;
		m_edgeElement.arrayNum(1) = 2 * m_edgeElement.nv_1 + 1;
		
		m_edgeElement.arrayNum(2) = 2 * m_edgeElement.nv_2 + 0;
		m_edgeElement.arrayNum(3) = 2 * m_edgeElement.nv_2 + 1;
		
		v_edgeElement.push_back(m_edgeElement);

		edgeNum = edgeNum + 1;
	}
	
}

void elasticPlate::computeBending()
{
	
	bendingNum = 0;
	v_bendingElement.clear();

	for (int i = 0; i < bending.size(); i++)
	{
		Vector3i bendingCurrent = bending[i];

		bendingElement m_bendingElement;

		m_bendingElement.nv_1 = bendingCurrent(0);
		m_bendingElement.nv_2 = bendingCurrent(1);
		m_bendingElement.nv_3 = bendingCurrent(2);

		m_bendingElement.x_1 = getVertex(m_bendingElement.nv_1);
		m_bendingElement.x_2 = getVertex(m_bendingElement.nv_2);
		m_bendingElement.x_3 = getVertex(m_bendingElement.nv_3);

		m_bendingElement.e_1 = m_bendingElement.x_2 - m_bendingElement.x_1;
		m_bendingElement.e_2 = m_bendingElement.x_3 - m_bendingElement.x_2;

		m_bendingElement.norm_1 =  m_bendingElement.e_1.norm();
		m_bendingElement.norm_2 =  m_bendingElement.e_2.norm();

		m_bendingElement.voroniLength = (m_bendingElement.norm_1 + m_bendingElement.norm_2) / 2;

		m_bendingElement.t_1 = m_bendingElement.e_1 / m_bendingElement.norm_1;
		m_bendingElement.t_2 = m_bendingElement.e_2 / m_bendingElement.norm_2;

		m_bendingElement.nBar = computeCurvature(m_bendingElement.x_1(0), m_bendingElement.x_1(1), m_bendingElement.x_2(0), m_bendingElement.x_2(1), m_bendingElement.x_3(0), m_bendingElement.x_3(1));

		m_bendingElement.arrayNum = VectorXi::Zero(6);

		m_bendingElement.arrayNum(0) = 2 * m_bendingElement.nv_1 + 0;
		m_bendingElement.arrayNum(1) = 2 * m_bendingElement.nv_1 + 1;
		
		m_bendingElement.arrayNum(2) = 2 * m_bendingElement.nv_2 + 0;
		m_bendingElement.arrayNum(3) = 2 * m_bendingElement.nv_2 + 1;
		
		m_bendingElement.arrayNum(4) = 2 * m_bendingElement.nv_3 + 0;
		m_bendingElement.arrayNum(5) = 2 * m_bendingElement.nv_3 + 1;

		if (i < curveNv - 2)
		{
			m_bendingElement.EI_local = EI;
		}
		else
		{
			m_bendingElement.EI_local = 1000 * EI;
		}

		v_bendingElement.push_back(m_bendingElement);

		bendingNum = bendingNum + 1;
	}
	
}



void elasticPlate::addElement()
{
	double minDistance = 1000.0;

	Vector2d xS = getVertex(0);
	int midIndex1 = 0;
	for (int i = curveNv; i < nv - 1; i++)
	{
		Vector2d xCurrent = ( getVertex(i) + getVertex(i+1) ) / 2;

		if ( (xS - xCurrent).norm() < minDistance )
		{
			minDistance = (xS - xCurrent).norm();
			midIndex1 = i;
		}
	}

	Vector2d xConnect = ( getVertex(midIndex1) + getVertex(midIndex1+1) ) / 2;

	Vector3i couplingLocal_1;
	couplingLocal_1(0) = midIndex1 + 2;
	couplingLocal_1(1) = midIndex1 + 1;
	couplingLocal_1(2) = 0;
	bending.push_back(couplingLocal_1);
	couplingLocal_1(0) = midIndex1 + 1;
	couplingLocal_1(1) = 0;
	couplingLocal_1(2) = 1;
	bending.push_back(couplingLocal_1);


	Vector2i edgeCurrent1;
	edgeCurrent1(0) = 0;
	edgeCurrent1(1) = midIndex1;
	edge.push_back(edgeCurrent1);
	edgeCurrent1(0) = 0;
	edgeCurrent1(1) = midIndex1 + 1;
	edge.push_back(edgeCurrent1);


	//cout << "xS " << xS.transpose() << endl;
	//cout << "xConnect " << xConnect.transpose() << endl;


	double offset = xConnect(1) - xS(1);

	for (int i = curveNv; i < nv; i++)
	{
		x(2 * i + 1) = x(2 * i + 1) - offset;
	}
	x0 = x;


	minDistance = 1000.0;

	Vector2d xE = getVertex(curveNv - 1);
	int midIndex2 = 0;
	for (int i = curveNv; i < nv - 1; i++)
	{
		Vector2d xCurrent = ( getVertex(i) + getVertex(i+1) ) / 2;

		if ( (xE - xCurrent).norm() < minDistance )
		{
			minDistance = (xE - xCurrent).norm();
			midIndex2 = i + 1;
		}
	}

	Vector2d xConnect2 = ( getVertex(midIndex2) + getVertex(midIndex2-1) ) / 2;

	Vector3i couplingLocal_2;
	couplingLocal_2(0) = midIndex2 - 2;
	couplingLocal_2(1) = midIndex2 - 1;
	couplingLocal_2(2) = curveNv - 1;
	bending.push_back(couplingLocal_2);
	couplingLocal_2(0) = midIndex2 - 1;
	couplingLocal_2(1) = curveNv - 1;
	couplingLocal_2(2) = curveNv - 2;
	bending.push_back(couplingLocal_2);


	Vector2i edgeCurrent2;
	edgeCurrent2(0) = curveNv - 1;
	edgeCurrent2(1) = midIndex2;
	edge.push_back(edgeCurrent2);
	edgeCurrent2(0) = curveNv - 1;
	edgeCurrent2(1) = midIndex2 - 1;
	edge.push_back(edgeCurrent2);

	//cout << "xE " << xE.transpose() << endl;
	//cout << "xConnect2 " << xConnect2.transpose() << endl;





	for (int i = edgeNum; i < edge.size(); i++)
	{
		Vector2i edgeCurrent = edge[i];

		edgeElement m_edgeElement;

		m_edgeElement.nv_1 = edgeCurrent(0);
		m_edgeElement.nv_2 = edgeCurrent(1);

		m_edgeElement.x_1 = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2 = getVertex(m_edgeElement.nv_2);

		m_edgeElement.refLength = (m_edgeElement.x_2- m_edgeElement.x_1).norm();
		m_edgeElement.edgeLength = m_edgeElement.refLength;

		//cout << m_edgeElement.nv_1 << " " << m_edgeElement.nv_2 << " " << m_edgeElement.refLength << endl;

		m_edgeElement.arrayNum = VectorXi::Zero(4);

		m_edgeElement.arrayNum(0) = 2 * m_edgeElement.nv_1 + 0;
		m_edgeElement.arrayNum(1) = 2 * m_edgeElement.nv_1 + 1;
		
		m_edgeElement.arrayNum(2) = 2 * m_edgeElement.nv_2 + 0;
		m_edgeElement.arrayNum(3) = 2 * m_edgeElement.nv_2 + 1;
		
		v_edgeElement.push_back(m_edgeElement);
	}

	//cout << edgeNum << " " << v_edgeElement.size() << endl;

	edgeNum = v_edgeElement.size();





	for (int i = bendingNum; i < bending.size(); i++)
	{
		Vector3i bendingCurrent = bending[i];

		bendingElement m_bendingElement;

		m_bendingElement.nv_1 = bendingCurrent(0);
		m_bendingElement.nv_2 = bendingCurrent(1);
		m_bendingElement.nv_3 = bendingCurrent(2);

		m_bendingElement.x_1 = getVertex(m_bendingElement.nv_1);
		m_bendingElement.x_2 = getVertex(m_bendingElement.nv_2);
		m_bendingElement.x_3 = getVertex(m_bendingElement.nv_3);

		m_bendingElement.e_1 = m_bendingElement.x_2 - m_bendingElement.x_1;
		m_bendingElement.e_2 = m_bendingElement.x_3 - m_bendingElement.x_2;

		m_bendingElement.norm_1 =  m_bendingElement.e_1.norm();
		m_bendingElement.norm_2 =  m_bendingElement.e_2.norm();

		m_bendingElement.voroniLength = (m_bendingElement.norm_1 + m_bendingElement.norm_2) / 2;

		m_bendingElement.t_1 = m_bendingElement.e_1 / m_bendingElement.norm_1;
		m_bendingElement.t_2 = m_bendingElement.e_2 / m_bendingElement.norm_2;

		m_bendingElement.nBar = computeCurvature(m_bendingElement.x_1(0), m_bendingElement.x_1(1), m_bendingElement.x_2(0), m_bendingElement.x_2(1), m_bendingElement.x_3(0), m_bendingElement.x_3(1));

		m_bendingElement.arrayNum = VectorXi::Zero(6);

		m_bendingElement.arrayNum(0) = 2 * m_bendingElement.nv_1 + 0;
		m_bendingElement.arrayNum(1) = 2 * m_bendingElement.nv_1 + 1;
		
		m_bendingElement.arrayNum(2) = 2 * m_bendingElement.nv_2 + 0;
		m_bendingElement.arrayNum(3) = 2 * m_bendingElement.nv_2 + 1;
		
		m_bendingElement.arrayNum(4) = 2 * m_bendingElement.nv_3 + 0;
		m_bendingElement.arrayNum(5) = 2 * m_bendingElement.nv_3 + 1;

		m_bendingElement.EI_local = 1000 * EI;

		v_bendingElement.push_back(m_bendingElement);

	}

	bendingNum = v_bendingElement.size();



}

double elasticPlate::computeCurvature(double xkm1, double ykm1, double xk, double yk, double xkp1, double ykp1)
{
	return 0.2e1 * tan(atan((-(xk - xkm1) * (ykp1 - yk) + (yk - ykm1) * (xkp1 - xk)) / ((xk - xkm1) * (xkp1 - xk) + (yk - ykm1) * (ykp1 - yk))) / 0.2e1);
}

void elasticPlate::clearMap()
{
	delete [] unconstrainedMap;
	delete [] fullToUnconsMap;

    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

void elasticPlate::reSetupMap()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;	
	
	// Setup the map from free dofs to all dof
	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];

	setupMap();
}

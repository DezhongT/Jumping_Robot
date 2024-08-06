#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include <fstream>

struct edgeElement
{
	int nv_1;
	int nv_2;

	Vector2d x_1;
	Vector2d x_2;

	double refLength;
	double edgeLength;

	VectorXi arrayNum;
};

struct bendingElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	VectorXi arrayNum;

	Vector2d x_1;
	Vector2d x_2;
	Vector2d x_3;

	Vector2d e_1;
	Vector2d e_2;

	double norm_1;
	double norm_2;

	Vector2d t_1;
	Vector2d t_2;

	double nBar;

	double voroniLength;

	double EI_local;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_radius, 
		double m_Possion, double m_dt, double m_deltaLength,
		 double m_l1, double m_l2, double m_compressRatio, double m_h1, double m_mass);
	~elasticPlate();

	double YoungM;
	double radius;
	double Possion;
	double dt;
	double density;

	double totalMass;

	Vector2d getVertex(int i);
	Vector2d getVertexStart(int i);
	Vector2d getVertexOld(int i);
	Vector2d getVelocity(int i);
	Vector2d getVelocityOld(int i);

	VectorXd x;
	VectorXd xSta;
	VectorXd x0;
	VectorXd u;

	std::vector<Vector2d> v_nodes;
    std::vector<Vector2i> edge;
    std::vector<Vector3i> bending;
    std::vector<int> constraint;

	std::vector<edgeElement> v_edgeElement;
	std::vector<bendingElement> v_bendingElement;

	double totalM;
	
	int temp;

	int nv;
	int edgeNum;
	int bendingNum;
	
	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector2d position, int k);

	void computeEdge();
	void computeBending();

	void updateEdgePair();
	void updateBendingPair();

	// boundary conditions
	int* isConstrained;
	int getIfConstrained(int k);
	int* unconstrainedMap;
	int* fullToUnconsMap;
	void setup();
	void setupMap();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
	void prepareForIteration();

	VectorXd massArray;
	void setupMass();
	VectorXi boundaryIndex;

	double EA;
	double EI;
	double crossSectionalArea;

	int curveNv;
	int frameNvS;
	int frameNvE;

	double computeCurvature(double xkm1, double ykm1, double xk, double yk, double xkp1, double ykp1);

	double l1;
    double l2;
    double compressRatio;
    double beaml;
    double h1;

    double deltaL;

    void clearMap();
    void reSetupMap();

    void addElement();

	private:
};

#endif

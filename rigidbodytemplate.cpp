#include "rigidbodytemplate.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <vector>
#include <Eigen/Geometry>       //cross product
#include <unordered_set>

using namespace std;
using namespace Eigen;

// LibIgl includes
#include <igl/doublearea.h>
#include <igl/centroid.h>

// RigidBodyTemplate Constructor
RigidBodyTemplate::RigidBodyTemplate(const MatrixX3d &verts, const MatrixX3i &faces) : V(verts), F(faces), area_(0)
{
    initialize();
}

// RigidBodyTemplate Destructor 
RigidBodyTemplate::~RigidBodyTemplate()
{    
	std::cout << "Rigid Body template memory resources being freed:\n";
}

// RigidBodyTemplate Initialize 
void RigidBodyTemplate::initialize()
{
    cout << "Executing RigidBodyTemplate Initialization code\n";
    // We already possess the list of boundary faces, of the rigid body [F]

    // 1) the volume of the body
	Eigen::VectorXd dblA;
	dblA.setZero();
	igl::doublearea(V,F,dblA);
	area_ = 0.5 * dblA.sum();

    // 2) the Center of Mass/Centroid of the body
    Eigen::Vector3d CoM(0,0,0);
	double area;
	igl::centroid(V,F,CoM,area);

    // translate the body, s.t. its center of mass lies at the origin (0,0,0)
	// #TODO :: is this translation of CoM coordinates correct? 2x-check here!
    oldCoM = CoM;
    for(int i = 0; i < V.rows(); ++i)
        V.row(i) -= CoM;
}

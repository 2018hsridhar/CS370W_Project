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

// source of confusion ... forces need to be scaled, according to the boundary length ( the boundary edge they apply too ). How to get this normal though? ... if it were vertex data, this would prove far far easier! 

// RigidBodyTemplate Constructor
RigidBodyTemplate::RigidBodyTemplate(const MatrixX3d &verts, const MatrixX3i &faces) : V(verts), F(faces), area_(0)
{
    inertiaTensor_.setZero();
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
	double vol;
	igl::centroid(V,F,CoM,vol);

    // translate the body, s.t. its center of mass lies at the origin (0,0,0)
    oldCoM = CoM;
    for(int i = 0; i < V.rows(); ++i)
        V.row(i) -= CoM;

    // 3) the intertia tensor (I do this for you here)

    // assumes you have computed F in step 1, and have translated the body to the origin
    computeInertiaTensor();
    cout << "DONE with inertia tensor" << endl;
}

// wait ... why does this compute normals too? Should I be worried about the normals computed for the inertia tensor? 
void RigidBodyTemplate::computeInertiaTensor()
{
    Vector3d quads(0,0,0);
    Vector3d mixed(0,0,0);
    for(int i=0; i<F.rows(); i++)
    {
        Vector3d pts[3];
        for(int j=0; j<3; j++)
        {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1]-pts[0]).cross(pts[2]-pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        Vector3d term(0,0,0);
        Vector3d mixterm(0,0,0);
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                for(int l=k; l<3; l++)
                {
                    for(int m=l; m<3; m++)
                        term[j] += pts[k][j]*pts[l][j]*pts[m][j];
                }
            }
            term[j] *= area*normal[j]/30.0;
        }
        double mix = 0;
        for(int j=0; j<3; j++)
        {
            mix += 6.0*pts[j][0]*pts[j][1]*pts[j][2];
            for(int k=0; k<3; k++)
            {
                mix += 2.0*pts[j][k]*pts[j][(k+1)%3]*pts[(j+1)%3][(k+2)%3];
                mix += 2.0*pts[j][k]*pts[j][(k+1)%3]*pts[(j+2)%3][(k+2)%3];
            }
            mix += pts[j][0]*pts[(j+1)%3][1]*pts[(j+2)%3][2];
            mix += pts[j][2]*pts[(j+1)%3][1]*pts[(j+2)%3][0];
        }
        for(int j=0; j<3; j++)
            mixterm[j] = mix*area*normal[j]/60.0;

        quads += term;
        mixed += mixterm;
    }

    inertiaTensor_ << quads[1]+quads[2], -mixed[2], -mixed[1],
            -mixed[2], quads[0]+quads[2], -mixed[0],
            -mixed[1], -mixed[0], quads[0]+quads[1];
}

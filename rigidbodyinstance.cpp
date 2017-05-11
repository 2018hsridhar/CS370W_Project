#include "rigidbodyinstance.h"
#include "vectormath.h"
#include "rigidbodytemplate.h"
#include <Eigen/Geometry>
#include <iostream>

using namespace Eigen;
using namespace std;

RigidBodyInstance::RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta) : c_0(c), c(c), theta_0(theta), theta(theta), rbtemplate_(rbtemplate)
{
    cvel.setZero();
    w.setZero();
}


#include "rigidbodyinstance.h"
#include "vectormath.h"
#include "rigidbodytemplate.h"
#include <Eigen/Geometry>
#include <iostream>

using namespace Eigen;
using namespace std;

RigidBodyInstance::RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta) : c(c), theta(theta), rbtemplate_(rbtemplate)
{
    cvel.setZero();
    w.setZero();
}


#ifndef RIGIDBODYINSTANCE_H
#define RIGIDBODYINSTANCE_H

#include <Eigen/Core>
#include <list>

class RigidBodyTemplate;

class RigidBodyInstance
{
public:
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta); 

	// NEED {c_0,theta_0} to capture changes :: {c_0,c} and {theta_0,theta}
    Eigen::Vector3d c_0;
    Eigen::Vector3d theta_0;

    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;

    const RigidBodyTemplate &getTemplate() {return rbtemplate_;}

private:
    const RigidBodyTemplate &rbtemplate_;
};

#endif // RIGIDBODYINSTANCE_H


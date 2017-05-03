#ifndef RIGIDBODYINSTANCE_H
#define RIGIDBODYINSTANCE_H

#include <Eigen/Core>
#include <list>

class RigidBodyTemplate;

class RigidBodyInstance
{
public:
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta); 

// I do need {c,theta} ... but aren't {cvel,w} always set to zero anyways? This kinda seems ... poitnless! 
    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;

    const RigidBodyTemplate &getTemplate() {return rbtemplate_;}

private:
    const RigidBodyTemplate &rbtemplate_;
};

#endif // RIGIDBODYINSTANCE_H

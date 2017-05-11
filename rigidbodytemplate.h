#ifndef RIGIDBODYTEMPLATE_H
#define RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>
#include <utility>
#include <functional>   //hash

class SignedDistanceField;

class RigidBodyTemplate
{
public:
    Eigen::Vector3d oldCoM;

    RigidBodyTemplate(const Eigen::MatrixX3d &verts, const Eigen::MatrixX3i &faces);
    ~RigidBodyTemplate();

    double getArea() const {return area_;}

    const Eigen::MatrixX3d &getVerts() const {return V;}
    const Eigen::MatrixX3i &getFaces() const {return F;}

private:
    RigidBodyTemplate(const RigidBodyTemplate &other);

    void initialize();

    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;

    double area_;

};

#endif // RIGIDBODYTEMPLATE_H

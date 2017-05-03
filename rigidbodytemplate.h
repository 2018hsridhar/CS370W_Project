#ifndef RIGIDBODYTEMPLATE_H
#define RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>
#include <utility>
#include <functional>   //hash

class SignedDistanceField;
// I also have the feeling that class SignedDistanceField, is basically pointless here? ... or perhaps not ( need distance to vertices, right? Any constraints for the Rigid Bodies needed too? )  

// set of SORTED faces, for the rigid body ( only given vertex-tet data?? huh?? ) well, i have both vertex and face data ... so a lot of this stuff, does NOT seem useful at all, for the template's set up? To focus on just using LibIgl calls in-place of this instead? Probably!



class RigidBodyTemplate
{
public:
    Eigen::Vector3d oldCoM;

    RigidBodyTemplate(const Eigen::MatrixX3d &verts, const Eigen::MatrixX3i &faces);
    ~RigidBodyTemplate();

    double getArea() const {return area_;}
    const Eigen::Matrix3d getInertiaTensor() const {return inertiaTensor_;}    

//    double getBoundingRadius() const {return radius_;}
    const Eigen::MatrixX3d &getVerts() const {return V;}
    const Eigen::MatrixX3i &getFaces() const {return F;}

private:
    RigidBodyTemplate(const RigidBodyTemplate &other);
    //RigidBodyTemplate &operator=(const RigidBodyTemplate &other);

    void initialize();

    void computeInertiaTensor();

    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;

    double area_;

   // double radius_;
	// calculate inertia tensor ... requires only face data! 
	// radius makes no sense, in our set up ( to set up a bounding circle? WHY ? )
	// wait ... do I need the inertia tensor? FOCUS on the matrices needed for F_ext to F_config! 
    Eigen::Matrix3d inertiaTensor_;
};

#endif // RIGIDBODYTEMPLATE_H

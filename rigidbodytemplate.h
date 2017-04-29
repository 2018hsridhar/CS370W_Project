#ifndef RIGIDBODYTEMPLATE_H
#define RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>
#include <utility>
#include <functional>   //hash

class SignedDistanceField;

class SortedFace
{
public:
    Eigen::Vector3d f;
    //int tIndex;
    int fIndex; //0 to 3; represents the first ijkl that was used in the set of 3

    //SortedFace(const Eigen::Vector3d& face, int tIndex, int fIndex);
    bool operator== (const SortedFace& other) const;
};

namespace std
{
    template <>
    struct hash<SortedFace>
    {
        //For our hash function
        size_t operator()(const SortedFace& sf) const
        {
            size_t seed = 0;
            for(int i = 0; i < 3; ++i)
            {
                double element = sf.f(i);
                seed ^= hash<double>()(element) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}

class RigidBodyTemplate
{
public:
    Eigen::Vector3d oldCoM;

 //   RigidBodyTemplate(const std::string &meshFilename, double scale);
    //RigidBodyTemplate(const Eigen::MatrixX3d &verts, const Eigen::MatrixX4i &tets);
    RigidBodyTemplate(const Eigen::MatrixX3d &verts);

    ~RigidBodyTemplate();

    double getVolume() const {return volume_;}
    const Eigen::Matrix3d getInertiaTensor() const {return inertiaTensor_;}    

    double getBoundingRadius() const {return radius_;}
    const Eigen::MatrixX3d &getVerts() const {return V;}
//    const Eigen::MatrixX4i &getTets() const {return T;}
    const Eigen::MatrixX3i &getFaces() const {return F;}

private:
    RigidBodyTemplate(const RigidBodyTemplate &other);
    //RigidBodyTemplate &operator=(const RigidBodyTemplate &other);

    void initialize();

    void computeInertiaTensor();

    Eigen::MatrixX3d V;
    //Eigen::MatrixX4i T;
    Eigen::MatrixX3i F;

    double volume_;
   // double radius_;
	// calculate inertia tensor ... requires only face data! 
	// radius makes no sense, in our set up! And neitehr does volume ( why does this matter ... for masses of the system? ) 
	// wait ... do I need the inertia tensor? FOCUS on the matrices needed for F_ext to F_config! 
    Eigen::Matrix3d inertiaTensor_;
};

#endif // RIGIDBODYTEMPLATE_H

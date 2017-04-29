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

/*
size_t Face_Hash::operator()(const SortedFace& sf) const
{
    return 0;
}
*/

bool SortedFace::operator== (const SortedFace& other) const
{
    const Eigen::Vector3d& l = this->f;
    const Eigen::Vector3d& r = other.f;
    return (l(0) == r(0)) && (l(1) == r(1)) && (l(2) == r(2));
}

SortedFace::SortedFace(const Eigen::Vector3d& face, int tIndex, int fIndex)
            : f(face), tIndex(tIndex), fIndex(fIndex)
{
    if(f(0) > f(1))
        std::swap(f(0), f(1));
    if(f(1) > f(2)) {
        std::swap(f(1), f(2));
        if(f(0) > f(1))
            std::swap(f(0), f(1));
    }
}

RigidBodyTemplate::RigidBodyTemplate(const std::string &meshFilename, double scale) : volume_(0), radius_(0)
{
    inertiaTensor_.setZero();

    ifstream ifs(meshFilename.c_str());
    vector<Vector3d> pts;
    vector<Vector4i> tets;
    while(true)
    {
        char c;
        ifs >> c;
        if(!ifs)
            break;
        if(c == 'v')
        {
            Vector3d pt;
            ifs >> pt[0] >> pt[1] >> pt[2];
            pts.push_back(pt);
        }
        if(c == 't')
        {
            Vector4i tet;
            ifs >> tet[0] >> tet[1] >> tet[2] >> tet[3];
            tets.push_back(tet);
        }
    }
    V.resize(pts.size(), 3);
    for(int i=0; i<(int)pts.size(); i++)
        V.row(i) = pts[i];
    T.resize(tets.size(), 4);
    for(int i=0; i<(int)tets.size(); i++)
    {
        T.row(i) = tets[i];
    }

    V *= scale;

    initialize();
}

RigidBodyTemplate::RigidBodyTemplate(const MatrixX3d &verts, const MatrixX4i &tets) : V(verts), T(tets), volume_(0)
{
    inertiaTensor_.setZero();
    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate()
{    
}

void RigidBodyTemplate::initialize()
{
    //cout << "Step 1" << endl;
    // need to compute the following things:
    // 1) the list of boundary faces of the rigid body, F
    typedef unordered_set<SortedFace/*, size_t (const SortedFace&), Face_Hash*/> MySet;
    MySet faceSet;
    for(int r = 0; r < T.rows(); ++r)
    {
        for(int fIndex = 0; fIndex < 4; ++fIndex)
        {
            //get our face
            Eigen::Vector3d f(T(r, (fIndex) % 4), T(r, (fIndex + 1) % 4), T(r, (fIndex + 2) % 4));
            SortedFace sf(f, r, fIndex);
            pair<MySet::iterator, bool> prevElement = faceSet.insert(sf);
            //if it was already in the set, remove it since it's not a boundary face
            if(!prevElement.second)
            {
//                std::cout << "old element" << endl;
                faceSet.erase(prevElement.first);
            }
            else
            {
//                std::cout << "NEW ELEMENT" << endl;
            }
        }
    }
    //resize F, & put in values
    F.resize(faceSet.size(), 3);
    int y = 0;
    for(auto it = faceSet.begin(); it != faceSet.end(); ++it)
    {
        const SortedFace& f = *it;
        int tIndex = f.tIndex;
        int fIndex = f.fIndex;
        Eigen::Vector3i trueF;
        //get outward-pointing face
        switch(fIndex)
        {
            case 0:
                trueF = Eigen::Vector3i(T(tIndex, 0), T(tIndex, 2), T(tIndex, 1));
                break;
            case 1:
                trueF = Eigen::Vector3i(T(tIndex, 1), T(tIndex, 2), T(tIndex, 3));
                break;
            case 2:
                trueF = Eigen::Vector3i(T(tIndex, 2), T(tIndex, 0), T(tIndex, 3));
                break;
            case 3:
                trueF = Eigen::Vector3i(T(tIndex, 3), T(tIndex, 0), T(tIndex, 1));
                break;
            default:
                //oops
                assert(false);
        }
        //add it to F
        for(int x = 0; x < 3; ++x)
            F(y,x) = trueF(x);
        //iterate
        ++y;
    }
    // 2) the volume of the body
    volume_ = 0;
    for(int i = 0; i < T.rows(); ++i)
    {
        const Eigen::Vector3d& pi = V.row(T(i,0));
        const Eigen::Vector3d& pj = V.row(T(i,1));
        const Eigen::Vector3d& pk = V.row(T(i,2));
        const Eigen::Vector3d& pl = V.row(T(i,3));
        volume_ += (pj - pi).cross(pk - pi).dot(pl - pi) / 6;
    }
    // 3) the center of mass of the body
    Eigen::Vector3d CoM(0,0,0);
    //for all faces, find xpart, ypart, and zpart
    for(int f = 0; f < F.rows(); ++f)
    {
        Eigen::Vector3i indices = F.row(f);
        Eigen::Vector3d vi = V.row(indices(0));
        Eigen::Vector3d vj = V.row(indices(1));
        Eigen::Vector3d vk = V.row(indices(2));
        //No need to normalize, since we'd just multiply by the length anyway
        Eigen::Vector3d n = (vj - vi).cross(vk - vi);
        //get components
        Eigen::Vector3d faceCoM;
        for(int q = 0; q < 3; ++q)
        {
            double i = vi(q);
            double j = vj(q);
            double k = vk(q);
            faceCoM(q) = (j * j + j * (k + i) + k * k + k * i + i * i) / 24;
        }
        CoM += faceCoM.cwiseProduct(n);
    }
    //correct by volume
    CoM /= volume_;
    //cout << "CoM: " << CoM(0) << " " << CoM(1) << " " << CoM(2) << endl;
    oldCoM = CoM;
    // then translate the body so that the center of mass is at (0,0,0)
    for(int i = 0; i < V.rows(); ++i)
        V.row(i) -= CoM;

    // 4) the radius of the smallest sphere enclosing the rigid body
    //only include vertices which exist in tets
    vector<bool> vertIsUsed(V.rows(), false);
    for(int i = 0; i < T.rows(); ++i)
    {
        const Eigen::Vector4i& t = T.row(i);
        for(int k = 0; k < 4; ++k)
            vertIsUsed[t(k)] = true;
    }
    for(int i = 0; i < V.rows(); ++i)
    {
        if(vertIsUsed[i])
        {
            double newR = V.row(i).norm();
            if(newR > radius_)
                radius_ = newR;
        }
    }
    //cout << "radius: " << radius_ << endl;
    // 5) the intertia tensor (I do this for you here)

    // assumes you have computed F in step 1, and have translated the body to the origin
    computeInertiaTensor();
    //cout << "done with inertia tensor" << endl;
}

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

#ifndef J_ALIGN_H
#define J_ALIGN_H

#include "tutorial_shared_path.h"
#include "vectormath.h"
#include <set>
#include <time.h> 

// LIBIGL includes [ some are and are not needed ]
#include <igl/random_dir.h> 
#include <igl/axis_angle_to_quat.h>
#include <igl/quat_to_mat.h>
#include <igl/doublearea.h>
#include <igl/per_vertex_normals.h>

// NEED for method headers 	
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
// J_ALIGN = Impulse based alignment of rigid bodies
using namespace Eigen; 

// #TODO :: ensure that none of these methods takes in a <mesh> struct as an input BTW !!!
namespace J_ALIGN 
{
	void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega);
	void solve_config_deltas(const Eigen::MatrixXd& vBoundary, const Eigen::MatrixXd bndExtForceMat, const RigidBodyInstance* body, Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega);
	void solveTransformation(const RigidBodyInstance* ScanBody, Eigen::Matrix4d& T);
	int mod(int a, int b);
	int detectFace(const Eigen::Vector2i& e_i, const Eigen::MatrixXi& F);
	void solveExternalForcesForBoundary(const std::vector<int>& boundaryLoop, const Eigen::MatrixXd& remeshV, const Eigen::MatrixXi& remeshF, Eigen::MatrixXd& boundaryV, Eigen::MatrixXd& extBndForceMat);
}
#endif

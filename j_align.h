#ifndef J_ALIGN_H
#define J_ALIGN_H

// includes for method header 
#include "rigidbodyinstance.h"
using namespace Eigen; 

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

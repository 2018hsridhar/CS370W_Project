#ifndef HELPERS_H
#define HELPERS_H

#include <Eigen/Core>
#include <Eigen/Geometry> 

using namespace Eigen; 

#include <iostream>
#include <vector>
#include <igl/cat.h>

namespace HELPER
{
	Eigen::MatrixXd convertToHomogenousForm(const Ref<const MatrixXd>& mat);
	Eigen::MatrixXd normalizeHomogenousMatrix (const Ref<const MatrixXd>& mat);
	void applyRigidTransformation(const Eigen::MatrixXd& V, const Eigen::MatrixXd& T, Eigen::MatrixXd& vOut);
	void constructNormalField(Eigen::MatrixXd& V, Eigen::MatrixXd& N, Eigen::MatrixXd& V_res, Eigen::MatrixXi& N_edges);
}

#endif 

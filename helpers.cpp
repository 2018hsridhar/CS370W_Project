#include "helpers.h"
#include <Eigen/Core>
#include <Eigen/Geometry> 

namespace HELPER 
{
	Eigen::MatrixXd convertToHomogenousForm(const Ref<const MatrixXd>& mat)
	{
		Eigen::MatrixXd homogenoizedMatrix;
		homogenoizedMatrix = mat.transpose().colwise().homogeneous().transpose(); 
		return homogenoizedMatrix;
	}

	Eigen::MatrixXd normalizeHomogenousMatrix (const Ref<const MatrixXd>& mat)
	{
		Eigen::MatrixXd normalizedMatrix;
		normalizedMatrix = mat.transpose().colwise().hnormalized().transpose(); 
		return normalizedMatrix;
	}

	void applyRigidTransformation(const Eigen::MatrixXd& V, const Eigen::MatrixXd& T, Eigen::MatrixXd& vOut)
	{
		Eigen::VectorXd zeroTranslate = Eigen::Vector4d (0,0,0,0);
		Eigen::MatrixXd p_i = HELPER::convertToHomogenousForm(V);
		Eigen::MatrixXd transV = (p_i * T).rowwise() + zeroTranslate.transpose();
		vOut = HELPER::normalizeHomogenousMatrix(transV);
	}
}

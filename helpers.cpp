#include <Eigen/Core>
#include <Eigen/Geometry> // for homogenous, and for hnormalized ... well, this is annoying!

//#include "glob_defs.h"
#include "helpers.h"

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
}

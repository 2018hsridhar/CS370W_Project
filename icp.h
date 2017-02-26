#ifndef ICP_H
#define ICP_H

#include "tutorial_shared_path.h"
using namespace Eigen;

namespace ICP
{
	void applyRigidICP(Eigen::MatrixXd& V_one, Eigen::MatrixXd& V_two, Eigen::MatrixXd& V_transformed);
}

#endif

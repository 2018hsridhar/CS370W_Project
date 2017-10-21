#ifndef BOOL_OPERS
#define BOOL_OPERS 

#include "tutorial_shared_path.h"
#include <Eigen/Core> 

namespace BOOL_OPERS 
{
	void applyUnitSphereRescaling(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_scaled);
}

#endif


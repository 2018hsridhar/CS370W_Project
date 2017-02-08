#ifndef HELPERS_H
#define HELPERS_H

#include <iostream> // lol forgot this! ... how does this keep messing me up!
#include <vector>
#include <Eigen/Core>

using namespace Eigen; 

namespace HELPER
{
	Eigen::MatrixXd convertToHomogenousForm(const Ref<const MatrixXd>& mat);
	Eigen::MatrixXd normalizeHomogenousMatrix (const Ref<const MatrixXd>& mat);
}

#endif // HELPERS_H




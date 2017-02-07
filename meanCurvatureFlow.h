#ifndef MCF_H
#define MCF_H

#include "tutorial_shared_path.h"
#include <Eigen/Core> 

namespace MCF
{
	bool meshHasBoundary(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F);
	void applyMeanCurvatureFlow(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F);
	void computeMeanCurvatureFlow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double timestep, Eigen::MatrixXd &Vc);
	void computeMeanCurvatureFlowWaterTight(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double timestep, Eigen::MatrixXd &Vc);
	void computeMeanCurvatureFlowBoundary(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, std::vector<bool> B, double timestep, Eigen::MatrixXd &Vc);

}

#endif


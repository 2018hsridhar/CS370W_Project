#ifndef BOOL_OPERS_H
#define BOOL_OPERS_H

#include "tutorial_shared_path.h"
#include <Eigen/Core> 

// #TODO :: convert to <namespace> pattern

namespace BOOL_OPERS 
{
	void applyUnitSphereRescaling(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_scaled, Eigen::MatrixXi& F_scaled);
	void throwMeshInSpace(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_p, Eigen::MatrixXi& F_p);
	void extractBoundary(Eigen::MatrixXd& VA, Eigen::MatrixXi& FA, Eigen::VectorXi& J, Eigen::MatrixXd& bV, Eigen::MatrixXi& bF);
	void applyBooleanOperations();
	void generateBoolOperMeshes(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V1, Eigen::MatrixXi& F1, Eigen::MatrixXd& V2, Eigen::MatrixXi& F2);
	void generateMeshFragment(Eigen::MatrixXd& V1, Eigen::MatrixXi& F1,
						Eigen::MatrixXd& V2, Eigen::MatrixXi& F2,
						Eigen::MatrixXd& fragV, Eigen::MatrixXi& fragF);
	void extractBoundary(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& J, Eigen::MatrixXd& BV, Eigen::MatrixXi& BF);

	// test logic < TO MOVE LATER >
	void testMeshThrow();
}

#endif


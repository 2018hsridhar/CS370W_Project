#ifndef BOOL_OPERS
#define BOOL_OPERS 

#include "tutorial_shared_path.h"
#include <Eigen/Core> 

// #TODO :: convert to <namespace> pattern

//namespace BOOL_OPERS 
//{
	void runIndexTestInterpRemesh();
	void applyUnitSphereRescaling(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V_scaled);
	void applyBooleanOperations();
	void generateBoolOperMeshes(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& V1, Eigen::MatrixXi& F1, Eigen::MatrixXd& V2, Eigen::MatrixXi& F2);
	void generateMeshFragment(Eigen::MatrixXd& V1, Eigen::MatrixXi& F1,
						Eigen::MatrixXd& V2, Eigen::MatrixXi& F2,
						Eigen::MatrixXd& fragV, Eigen::MatrixXi& fragF);
//}

#endif


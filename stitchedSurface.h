#ifndef INTERP_H 
#define INTERP_H 

#include "path.h"
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/cat.h> 
#include <igl/point_mesh_squared_distance.h>
#include <igl/slice.h>
#include <igl/boundary_loop.h>
using namespace Eigen; 

// #TODO :: how to fix this "using namespace" annoying issue? 

// #TODO  :: store a notion of the original markers? unsure atm
namespace STITCHED_SURF 
{
	int mod(int a, int b);
	//void retIndicesBoundaryVertices(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, std::vector<int>& bndIndices);
	//void generateBoundaryVertices(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd& V_boundary);
	void generateStitchedSurface(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1, 
								const Eigen::MatrixXd &V2, const Eigen::MatrixXi &F2, 
								Eigen::MatrixXd& vOff, Eigen::MatrixXi& fOff);
	bool edgesAreEqual(const Ref<const Eigen::Vector2i>& e1, const Ref<Eigen::Vector2i>& e2);
	void constructSeedEdgeIndices(Eigen::VectorXi& bndIndexesSrc, Eigen::MatrixXd& bndVertsSrc,
							Eigen::VectorXi& bndIndexesDest, Eigen::MatrixXd& bndVertsDest,	
							Eigen::Vector2i& seedEdge);
}
#endif

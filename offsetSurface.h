#ifndef OFFSET_H 
#define OFFSET_H 

#include "tutorial_shared_path.h"
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/cat.h> 
#include <igl/point_mesh_squared_distance.h>
#include <igl/slice.h>
#include <igl/boundary_loop.h>
using namespace Eigen; // how to fix this "using namespace" annoying issue? 

// #TODO  :: store a notion of the original markers? unsure atm
namespace OFFSET
{
	void generateOffsetSurface(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1, 
								const Eigen::MatrixXd &V2, const Eigen::MatrixXi &F2, 
								Eigen::MatrixXd& vOff, Eigen::MatrixXi& fOff);
	bool edgesAreEqual(const Ref<const Eigen::Vector2i>& e1, const Ref<Eigen::Vector2i>& e2);
	void constructSeedEdge(Eigen::VectorXi& bndIndexesV1, Eigen::MatrixXd& bndVertsV1,Eigen::VectorXi& bndVertsV2,Eigen::Vector2i& seedEdge);



//void fillBoundaryIndexes(std::vector<bool> verticesBoundaryStatus, int numBoundaryVertices,int* indexesArray);
//template <typename T>
//void printcoll (T const& coll);
//Eigen::MatrixXd retrieveBoundaryVerticesInScan1(int bndryVerticesIdxs[], int numBndryVertices);
//Eigen::MatrixXd retrieveBoundaryVerticesInScan2(int bndryVerticesIdxs[], int numBndryVertices);
//std::vector<int> findAdjBndryVertsInScan1(int vertIdx);
//std::vector<int> findAdjBndryVertsInScan2(int vertIdx);
//int findClosestAdjBndryNodeInScan1(std::vector<int> bndryVertices, int myVertexIndex, int prev_v1);
//int findClosestAdjBndryNodeInScan2(std::vector<int> bndryVertices, int myVertexIndex, int prev_v2);
}

#endif

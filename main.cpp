﻿// PURPOSE :: a greedy surface reconstruction algorithm, from 2 range images, based of closest distances between boundary vertices

#include "tutorial_shared_path.h"
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/cat.h> 
#include <igl/point_mesh_squared_distance.h>
#include <igl/slice.h>
#include <igl/boundary_loop.h>

using namespace Eigen;  
using namespace igl;

///////// FORWARD METHOD DECLARATIONS /////////////////////
//void fillBoundaryIndexes(std::vector<bool> verticesBoundaryStatus, int numBoundaryVertices,int* indexesArray);
//template <typename T>
//void printcoll (T const& coll);
//Eigen::MatrixXd retrieveBoundaryVerticesInScan1(int bndryVerticesIdxs[], int numBndryVertices);
//Eigen::MatrixXd retrieveBoundaryVerticesInScan2(int bndryVerticesIdxs[], int numBndryVertices);
//std::vector<int> findAdjBndryVertsInScan1(int vertIdx);
//std::vector<int> findAdjBndryVertsInScan2(int vertIdx);
//int findClosestAdjBndryNodeInScan1(std::vector<int> bndryVertices, int myVertexIndex, int prev_v1);
//int findClosestAdjBndryNodeInScan2(std::vector<int> bndryVertices, int myVertexIndex, int prev_v2);
bool edgesAreEqual(const Ref<const Eigen::Vector2i>& e1, const Ref<Eigen::Vector2i>& e2);

///////// PREALLOCATION ///////////////////
struct Mesh
{
  Eigen::MatrixXd V; 
  Eigen::MatrixXi F;
  Eigen::MatrixXi E; 
} scan1,scan2,scans,scene,interpolatedSurface;

std::vector<std::vector<int>> Adjacency_Scan1;
std::vector<std::vector<int>> Adjacency_Scan2;

std::vector<bool> boundaryVerticesStatus_scan1;
std::vector<bool> boundaryVerticesStatus_scan2;

Eigen::MatrixXd boundaryVertices_scan1;
Eigen::MatrixXd boundaryVertices_scan2;

int main(int argc, char *argv[])
{
	std::cout << "executing offset surface generation\n" << std::endl;

  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy.off",scan1.V,scan1.F)) {
    std::cout << "Failed to load partial scan one." << std::endl;
  } 
  if(!readOFF(TUTORIAL_SHARED_PATH "/planexy2.off",scan2.V,scan2.F)) {
    std::cout << "Failed to load partial scan two." << std::endl;
  }

// PREALLOC 
std::vector<Eigen::Vector2i> allEdges;
std::vector<int> newTriangleFaces;

///////// SOLVE BOUNDARY VERTICES ( GET CYCLICAL ORDERING TOO ) ///////
// note that the set of boundary vertex INDICES will be ordered, in a cyclical manner! 
// ... the real issue, is determing which cyclic order to work with ( might need to reverse bndIndexesScan2)
Eigen::VectorXi bndIndexesScan1;
igl::boundary_loop(scan1.F,bndIndexesScan1); 
Eigen::MatrixXd bndVertsScan1;
igl::slice(scan1.V,bndIndexesScan1,1,bndVertsScan1); // #TODO :: refactor?? 1 is for a row option, I believe! 
int numBoundaryVerticesScan1 = bndIndexesScan1.rows();

Eigen::VectorXi bndIndexesScan2;
igl::boundary_loop(scan2.F,bndIndexesScan2); 
Eigen::VectorXi revCycle = bndIndexesScan2.reverse();
Eigen::MatrixXd bndVertsScan2;
igl::slice(scan2.V,revCycle,1,bndVertsScan2);
int numBoundaryVerticesScan2 = bndIndexesScan2.rows();

/*
  [1] solve for a seed edge :: choose a rand point in scan_1, find closest point in scan_2
  [2] keep alternating edge solving , and use the adjacency lists  
  [3] end once you have the original edge data ! Add this last face 
 */

///////////////////////////////////////////////////////////////////////////////////////////////
// [1] solve for a seed edge :: choose a rand point in scan_1, find closest point in scan_2 ///
///////////////////////////////////////////////////////////////////////////////////////////////
/// #TODO :: name this better!

  std::cout << "SETTING up a seed edge.\n ";
  int scan1SeedIdx = bndIndexesScan1(0); // i, @ this point, itself, is 0 
  Eigen::MatrixXd scan1SeedPoint = bndVertsScan1.row(0);

  Eigen::MatrixXd closestPointToSeedInScan2;
  Eigen::VectorXd smallestSquaredDists;
  Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(bndVertsScan2.rows(), 0, bndVertsScan2.rows() - 1);
  Eigen::VectorXi smallestDistIndxs;
  igl::point_mesh_squared_distance(scan1SeedPoint,bndVertsScan2,
                                    Ele,
									smallestSquaredDists,smallestDistIndxs,
									closestPointToSeedInScan2);

  //std::cout << "Index = " << smallestDistIndxs(0) << std::endl;  
  int scan2ClosestPointToSeedIndex = bndIndexesScan2(smallestDistIndxs(0));
  Eigen::Vector2i seedEdge = Eigen::Vector2i( scan1SeedIdx, scan2ClosestPointToSeedIndex + scan1.V.rows());
 
  allEdges.push_back(seedEdge);
  Eigen::Vector2i newEdge;

  //std::cout << seedEdge.transpose() << std::endl;

// [2] keep alternating edge solving , and use the adjacency lists  
	int i = 0; 
	int j = smallestDistIndxs(0);

    std::cout << "PROGRESSING over the Greedy Zippering Surface Reconstruction Algorithm.\n ";
    do
    {
	/*
	 NOTE where access to information lies ... index into boundary vertices, vs bndry vertices itself, vs index!
		- int bndryVertexIdx = bndIndexesScan1(i); // i, @ this point, itself, is 0 
		- Eigen::MatrixXd bndryVertexNode = bndVertsScan1.row(i);
	*/

		// #TODO :: fix bug in i-j update
   		//std::cout << " indices i and j are \t\t " << i << "\t\t" << j << std::endl;
		int p_i = bndIndexesScan1(i);
        //std::cout << "can access p_i" << std::endl;
		int q_j = bndIndexesScan2(j);
        //std::cout << "can access q_j" << std::endl; // ERR :: cannot access this! 

        int p_i_plus_1 = bndIndexesScan1((i + 1) % numBoundaryVerticesScan1);
        int q_j_plus_1 = bndIndexesScan2((j + 1) % numBoundaryVerticesScan2); 

    	//std::cout << " bndry vertices p_i, q_j are \t\t " << p_i << "\t\t" << q_j << std::endl;
    	//std::cout << " the next bndry vertices p_i_plus_1, q_j_plus_1, are \t\t " << p_i_plus_1 << "\t\t" << q_j_plus_1 << std::endl;
 		//  assess d(b_p_i, b_p_i_plus_1) <= d(b_q_j, b_q_j_plus_1)
	    Eigen::VectorXd b_p_i = scan1.V.row(p_i);
		Eigen::VectorXd b_q_j = scan2.V.row(q_j);

		Eigen::VectorXd b_p_i_plus_1 = scan1.V.row(p_i_plus_1);
		Eigen::VectorXd b_q_j_plus_1 = scan2.V.row(q_j_plus_1);
     
        //std::cout << "Can Access border vertex data\n";

		double scan1Distance =  (b_p_i - b_p_i_plus_1).squaredNorm();
		double scan2Distance =  (b_q_j - b_q_j_plus_1).squaredNorm();
		//std::cout << scan1Distance << '\t' << scan2Distance << std::endl;
		bool isItEdgeInScanOne = (scan1Distance < scan2Distance);
		// construct new edge
		Eigen::Vector3i newFace; 			// #TODO :: ensure that this is correct! 
		if (isItEdgeInScanOne) {
			newEdge = Eigen::Vector2i(p_i_plus_1, (q_j + scan1.V.rows()));
			newFace = Eigen::Vector3i(p_i,(q_j + scan1.V.rows()),p_i_plus_1);
			i = (i + 1) % numBoundaryVerticesScan1;
		}
		else {  
			newEdge = Eigen::Vector2i(p_i, (q_j_plus_1 + scan1.V.rows()));
			newFace = Eigen::Vector3i(p_i,(q_j + scan1.V.rows()),(q_j_plus_1+ scan1.V.rows()));
			j = (j + 1) % numBoundaryVerticesScan2;
		}

		// now that new edge is added, continue on with the algorithm ! 
		allEdges.push_back(newEdge); 
		newTriangleFaces.push_back(newFace(0));
		newTriangleFaces.push_back(newFace(1));
		newTriangleFaces.push_back(newFace(2));
    } while(!edgesAreEqual(newEdge,seedEdge));

	//////////////////////////////////////////////////////////////////////////
	// [3] end once you have the original edge data ! Add this last face  ///
	//////////////////////////////////////////////////////////////////////////
	// CONVERT the set of ( 3 * faces ) integers, of vertex indices, to a matrix ( for faces data )

  std::cout << "Constructing interpolating surface vertex and face data.\n";
  int numOfFaces = newTriangleFaces.size() / 3;
  Eigen::MatrixXi faces = Eigen::Map<Eigen::MatrixXi,RowMajor> (&newTriangleFaces[0],3,numOfFaces); // #TODO :: check if this is correct!
  interpolatedSurface.F = faces.transpose();
  std::cout << interpolatedSurface.F << std::endl;

  // CREATE ONE HUGE MESH containing the two partial scans and interpolated surface
  igl::cat(1,scan1.V,scan2.V,scene.V);
  igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), scans.F);
  igl::cat(1,scans.F, interpolatedSurface.F, scene.F);

/*
  // SOLVE FOR A TETRAHEDRALIZATION, based on (scene.V, scene.F)
	// triangulate makes sense ONLY for 2-dimensional data ( i.e. triangulate a polygon, or a boundary ) ... okay ... I need to study this differenec, for sure { when to triangulate, versus tetrahedralize, a mesh/surface } ... since I care about the boundary, and since the interpolating surface, between the boundary, doesn't need any more vertices insside ... it's a triangulation !


	// seems liike it shoudl work for an interp surface ( nothing inside ) ... but vertex is 3-d data... hence, kinda screwed 
	// but generating vertices inside ... that is a job best reserved for TETGEN!

	// do I care ONLY about the convex hull, or also the "interior" of a given surface?? not sure!
	// ... I guess a tetrahedralized interior! 
	Eigen::MatrixXd TV;
	Eigen::MatrixXi TT;
	Eigen::MatrixXi TF;
	/// ... something seems ammiss here /// 
	igl::tetrahedralize(scene.V,scene.F,"pq1.414a0.01", TV,TT,TF);	
*/
	

  /***********************************************************/ 
  // SETUP LibIgl Viewer 
  /***********************************************************/ 
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(scene.V, scene.F); 
  //viewer.data.set_mesh(TV,TF);
  viewer.launch();
 
}

bool edgesAreEqual(const Ref<const Eigen::Vector2i>& e1, const Ref<Eigen::Vector2i>& e2)
{
	bool edgesAreEqual = false;
	double edgeNormDist =  (e1- e2).norm();
	edgesAreEqual = (edgeNormDist == 0);
    return edgesAreEqual;
}


﻿// PURPOSE :: a greedy surface reconstruction algorithm, from 2 range images, based of closest distances between boundary vertices

#include "offsetSurface.h"
#include "remesh.h"
#include "glob_defs.h"
using namespace Eigen;  
using namespace igl;

// #TODO :: in glob_defs, have variables you can set ( for say, edge lengths in CGAl code, or files you want to test ... wait, for files you want to tet, that's a good @(#)$u idea ) . Your files should or should not be different, based on what you do and do not want to test. A.K.A. make it configurable!
// #TODO :: if time available, convert main() to a testing method

int main(int argc, char *argv[])
{
	// LOAD IN MESHES 
	std::cout << "Executing (TEST) for offset surface generation." << std::endl;
	struct Mesh
	{
		Eigen::MatrixXd V; 
		Eigen::MatrixXi F;
	} scan1,scan2,scene, remeshed;

	if(!readOFF(GLOBAL::offsetGenScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::offsetGenScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	OFFSET::generateOffsetSurface(scan1.V,scan1.F,scan2.V,scan2.F, scene.V,scene.F);
	//REMESH::remeshSurface(scene.V,scene.F,remeshed.V,remeshed.F);

	//  SETUP LibIgl Viewer 
	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(scene.V, scene.F); 
	//viewer.data.set_mesh(remeshed.V, remeshed.F); 
	viewer.launch();
}

namespace OFFSET
{
	void generateOffsetSurface(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1, 
								const Eigen::MatrixXd &V2, const Eigen::MatrixXi &F2, 
								Eigen::MatrixXd& vOff, Eigen::MatrixXi& fOff)
	{
		///////// PREALLOCATION ///////////////////
		// #TODO :: check if this is used later :-P 
		struct Mesh
		{
			Eigen::MatrixXd V; 
			Eigen::MatrixXi F;
		} scans,interpolatedSurface;

		std::vector<std::vector<int>> Adjacency_Scan1;
		std::vector<std::vector<int>> Adjacency_Scan2;

		std::vector<bool> boundaryVerticesStatus_scan1;
		std::vector<bool> boundaryVerticesStatus_scan2;

		Eigen::MatrixXd boundaryVertices_scan1;
		Eigen::MatrixXd boundaryVertices_scan2;

		std::vector<Eigen::Vector2i> allEdges;
		std::vector<int> newTriangleFaces;

		///////// SOLVE BOUNDARY VERTICES ( GET CYCLICAL ORDERING TOO ) ///////
		/*
 		 * note that the set of boundary vertex INDICES will be ordered, in a cyclical manner! 
		 * #TODO :: DETERMINE optional cycle order ( might need to reverse bndIndexesScan2)
		 */
		Eigen::VectorXi bndIndexesScan1;
		igl::boundary_loop(F1,bndIndexesScan1); 
		Eigen::MatrixXd bndVertsScan1;
		igl::slice(V1,bndIndexesScan1,1,bndVertsScan1); 
		int numBoundaryVerticesScan1 = bndIndexesScan1.rows();

		Eigen::VectorXi bndIndexesScan2;
		igl::boundary_loop(F2,bndIndexesScan2); 
		Eigen::MatrixXd bndVertsScan2;
		Eigen::VectorXi revCycle = bndIndexesScan2.reverse();
		igl::slice(V1,bndIndexesScan2,1,bndVertsScan2); 
	//	igl::slice(V2,revCycle,1,bndVertsScan2);
		int numBoundaryVerticesScan2 = bndIndexesScan2.rows();

		/*
	     * Algorithm @ a high level
	     * [1] solve for a seed edge :: choose a rand point in scan_1, find closest point in scan_2
		 * [2] keep alternating edge solving , and use the adjacency lists  
		 * [3] end once you have the original edge data ! Add this last face 
		 */

		/* 
		 * [1] solve for a seed edge :: 
		 * choose a rand point in scan_1, find closest point in scan_2 ///
		 */

		std::cout << "[1] SOLVE for the seed edge.\n";
	
		Eigen::Vector2i seedEdgeIndices;
		OFFSET::constructSeedEdgeIndices( bndIndexesScan1, bndVertsScan1, 
							bndIndexesScan2, bndVertsScan2, seedEdgeIndices); 

		Eigen:Vector2i seedEdge;
		int i = seedEdgeIndices(0); 
		int j = seedEdgeIndices(0);
		seedEdge(0) = bndIndexesScan1(i);
		seedEdge(1) =  bndIndexesScan2(j) + V1.rows();
		allEdges.push_back(seedEdge);
		Eigen::Vector2i newEdge;

		// [2] keep alternating edge solving, using the adj lists
		std::cout << "PROGRESSING over Greedy Zippering Surface Reconstruction Algorithm.\n ";
		// why is it that one of these vertices seems to always be fixed? that's weird
		int iter = 0; 
		do
		{
			int p_i = bndIndexesScan1(i);
			int q_j = bndIndexesScan2(j);

			int p_i_plus_1 = bndIndexesScan1((i + 1) % numBoundaryVerticesScan1);
			int q_j_plus_1 = bndIndexesScan2((j + 1) % numBoundaryVerticesScan2); 

			//  assess d(b_p_i, b_p_i_plus_1) <= d(b_q_j, b_q_j_plus_1)
			Eigen::VectorXd b_p_i = V1.row(p_i);
			Eigen::VectorXd b_q_j = V2.row(q_j);

			Eigen::VectorXd b_p_i_plus_1 = V1.row(p_i_plus_1);
			Eigen::VectorXd b_q_j_plus_1 = V2.row(q_j_plus_1);
		
			double scan1Distance =  (b_p_i - b_p_i_plus_1).squaredNorm();
			double scan2Distance =  (b_q_j - b_q_j_plus_1).squaredNorm();
			bool isItEdgeInScanOne = (scan1Distance < scan2Distance);

			// construct new edge
			Eigen::Vector3i newFace; 			
			if (isItEdgeInScanOne) {
				newEdge = Eigen::Vector2i(p_i_plus_1, (q_j + V1.rows()));
				newFace = Eigen::Vector3i(p_i,(q_j + V1.rows()),p_i_plus_1);
				i = (i + 1) % numBoundaryVerticesScan1;
			}
			else {  
				// always seesm to execute for weird camel head case?? WTF??
				newEdge = Eigen::Vector2i(p_i, (q_j_plus_1 + V1.rows()));
				//newFace = Eigen::Vector3i(p_i,(q_j + V1.rows()),(q_j_plus_1+ V1.rows()));
				newFace = Eigen::Vector3i(p_i,(q_j_plus_1 + V1.rows()),(q_j + V1.rows()));
				j = (j + 1) % numBoundaryVerticesScan2;
			}

			// push new edge, progress the algorithm
			allEdges.push_back(newEdge); 
			newTriangleFaces.push_back(newFace(0));
			newTriangleFaces.push_back(newFace(1));
			newTriangleFaces.push_back(newFace(2));
			//std::cout << "new edge = " << newEdge << std::endl;
			iter++;
		} while(iter < 10000 && !OFFSET::edgesAreEqual(newEdge,seedEdge));
		//} while(iter < 1000);
		//} while(!OFFSET::edgesAreEqual(newEdge,seedEdge));

		//////////////////////////////////////////////////////////////////////////
		// [3] end once you have the original edge data ! Add this last face  ///
		//////////////////////////////////////////////////////////////////////////
		// CONVERT set of ( 3 * faces ) integers, of vertex indices, to a matrix ( for faces data )

		std::cout << "Constructing interpolating surface vertex and face data.\n";
		int numOfFaces = newTriangleFaces.size() / 3;
		Eigen::MatrixXi faces = Eigen::Map<Eigen::MatrixXi,RowMajor> (&newTriangleFaces[0],3,numOfFaces); 
		interpolatedSurface.F = faces.transpose();

		// CREATE one huge interpolating mesh that contians 
		// both the two partial scans, and the interpolating surface
		igl::cat(1,V1,V2,vOff);
		igl::cat(1,F1, MatrixXi(F2.array() + V1.rows()), scans.F);
		igl::cat(1,scans.F, interpolatedSurface.F, fOff);
	}

	/*
     * method does not account for offset by vertices in a scan.
     * It assumes that for both scans, the boundary indices 
     * commence from 0 to scan length - 1
     */
	void constructSeedEdgeIndices(Eigen::VectorXi& bndIndexesSrc, Eigen::MatrixXd& bndVertsSrc,
							Eigen::VectorXi& bndIndexesDest, Eigen::MatrixXd& bndVertsDest,	
							Eigen::Vector2i& seedEdge)
	{
        int scan1SeedIdx = bndIndexesSrc(0); 
		Eigen::MatrixXd scan1SeedPoint = bndVertsSrc.row(0);

		Eigen::MatrixXd closestPointToSeedInScan2;
		Eigen::VectorXd smallestSquaredDists;
		Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(bndVertsDest.rows(), 0, bndVertsDest.rows() - 1);
		Eigen::VectorXi smallestDistIndxs;
		igl::point_mesh_squared_distance(scan1SeedPoint,bndVertsDest,
										Ele,
										smallestSquaredDists,smallestDistIndxs,
										closestPointToSeedInScan2);

		seedEdge = Eigen::Vector2i( 0, smallestDistIndxs(0));
	}	

	bool edgesAreEqual(const Ref<const Eigen::Vector2i>& e1, const Ref<Eigen::Vector2i>& e2)
	{
		bool edgesAreEqual = false;
		double edgeNormDist =  (e1- e2).norm();
		edgesAreEqual = (edgeNormDist == 0);
		return edgesAreEqual;
	}
}

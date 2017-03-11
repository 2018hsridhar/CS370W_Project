﻿// PURPOSE :: a greedy interpolating surface construction algorithm, from 2 range images, based of sticthing the best edge in a 4-point processing method. #TODO :: write up with a better description!

#include "interpSurface.h"
#include "remesh.h"
#include "glob_defs.h"
using namespace Eigen;  
using namespace igl;
using namespace std; // take out, but for now, useful in debug statements

//#include "igl/orient_outward.h"
// use combine.h --- for mesh concatenation?
// use cut_mesh.h --- for cutting? anyway to just cut along an axis ... really!! ... nope, need the set of edges to cut 
// use orient_outward.h --- to orient mesh outwards?? seems useful?? ... possible, might fail
// use slice methods too - only for matrices, or for tetrahedral meshes
// easiest -- rotate the plane, by a 180-deg rotation, about its y-axis ( rip from generate_data.cpp). This will prove easiest!

namespace INTERP_SURF 
{

	int mod(int a, int b)
	{
		return (a%b+b)%b;
	}

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
		std::vector<int> newTriangleFaces;

		/*
		 * *** SOLVE BOUNDARY VERTICES ( GET CYCLICAL ORDERING TOO ) ***
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
		igl::slice(V2,bndIndexesScan2,1,bndVertsScan2);  
		int numBoundaryVerticesScan2 = bndIndexesScan2.rows();

// trying to walk in wrong way ... orientation is random in boundry cycle, or how output is used
// they're walking the wrong ways!

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
		INTERP_SURF::constructSeedEdgeIndices( bndIndexesScan1, bndVertsScan1, 
							bndIndexesScan2, bndVertsScan2, seedEdgeIndices); 

		Eigen:Vector2i seedEdge;
		int i = seedEdgeIndices(0); 
		int j = seedEdgeIndices(1);
		
		seedEdge(0) = bndIndexesScan1(i);
		seedEdge(1) =  bndIndexesScan2(j) + V1.rows();
		Eigen::Vector2i newEdge;

		// [2] keep alternating edge solving, using the adj lists
	// NOTE :: you DO NOT walk in the same direction, for these cases!
		std::cout << "PROGRESSING over Greedy Zippering Surface Reconstruction Algorithm.\n ";
		int iter = 0; 
		bool edgesToScan1 = false;
		bool edgesToScan2 = false; 
		do
		{
			// [1] GATHER INDEX AND VERTEX DATA 
			// set up {p_i,p_i+1, q_i,q_i+1}
			int p_i = bndIndexesScan1(i);
			int q_j = bndIndexesScan2(j);

			int p_i_plus_1 = bndIndexesScan1((i + 1) % numBoundaryVerticesScan1);
			int q_j_plus_1 = bndIndexesScan2(INTERP_SURF::mod(j-1,numBoundaryVerticesScan2));

			//  assess d(b_p_i, b_p_i_plus_1) <= d(b_q_j, b_q_j_plus_1)
			Eigen::VectorXd b_p_i = V1.row(p_i);
			Eigen::VectorXd b_q_j = V2.row(q_j);

			Eigen::VectorXd b_p_i_plus_1 = V1.row(p_i_plus_1);
			Eigen::VectorXd b_q_j_plus_1 = V2.row(q_j_plus_1);
	
			// adding the edges, as usual
			double scan1Distance =  (b_q_j - b_p_i_plus_1).squaredNorm(); // from scan 1
			double scan2Distance =  (b_p_i - b_q_j_plus_1).squaredNorm(); // scan 2
			bool isItEdgeInScanOne = (scan1Distance < scan2Distance);
			Eigen::Vector3i newFace; 			

			if (isItEdgeInScanOne) 
			{
				newEdge = Eigen::Vector2i(p_i_plus_1, (q_j + V1.rows()));
				newFace = Eigen::Vector3i(q_j + V1.rows(), p_i_plus_1,p_i);
				i = (i + 1) % numBoundaryVerticesScan1;

				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));
		
				// check here, if you are back @ seed edge! must be here!
				if(p_i_plus_1 == seedEdge(0))
				{
					edgesToScan2 = true;
					std::cout << "Met same p_i seed edge point" << std::endl;
					break;
				}

				//	std::cout << "in scan 1 " << std::endl;
			}
			else 
			{  
				newEdge = Eigen::Vector2i(p_i, (q_j_plus_1 + V1.rows()));
				newFace = Eigen::Vector3i(p_i,(q_j + V1.rows()),(q_j_plus_1 + V1.rows()));
				j = INTERP_SURF::mod(j - 1, numBoundaryVerticesScan2);
			//	std::cout << "in scan 2 " << std::endl;

				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));

				if (q_j_plus_1 + V1.rows() == seedEdge(1)) 
				{
					edgesToScan1 = true;
					std::cout << "Met same q_j seed edge point" << std::endl;
					break;
				}
			}
	
			iter++;
			std::cout << iter << std::endl;

		} 
		while(!INTERP_SURF::edgesAreEqual(newEdge,seedEdge));
		//while(iter <= 0);
		//while(iter <= 200);

		// [2] CLEAN UP PHASE - if the new Edge happens to be adjacent to a seed vertex
		iter = 0;
		// do not push back 2 repeated triangles --- you have to look in the mesh, to detect this :-P!
		if(edgesToScan2)
		{
			std::cout << "Cleaning up to Scan 2" << '\n';
			assert(bndIndexesScan1(i) == seedEdge(0));
			int p_i = bndIndexesScan1(i);

			// #TODO :: check if do-while, mustbe converted to a while loop
			// I am basically overwriting my first triangle here ... that is the issue!
			do
			{
				int q_j = bndIndexesScan2(j);
				int q_j_plus_1 = bndIndexesScan2(INTERP_SURF::mod(j-1,numBoundaryVerticesScan2));

				newEdge = Eigen::Vector2i(p_i,(q_j_plus_1 + V1.rows()));
				Eigen::Vector3i newFace = Eigen::Vector3i(p_i,(q_j + V1.rows()),(q_j_plus_1 + V1.rows()));
				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));
				j = INTERP_SURF::mod(j - 1, numBoundaryVerticesScan2);
				cout << "iter = " << iter << endl;
				iter++;
			} while(!INTERP_SURF::edgesAreEqual(newEdge,seedEdge));
		}
		if (edgesToScan1)
		{
			std::cout << "Cleaning up to Scan 1" << '\n';
			assert(bndIndexesScan2(j) == seedEdge(1));
			int q_j = bndIndexesScan2(j);

			do
			{
				int p_i = bndIndexesScan1(i);
				int p_i_plus_1 = bndIndexesScan1((i+1) % numBoundaryVerticesScan1);

				newEdge = Eigen::Vector2i(p_i_plus_1, (q_j + V1.rows()));
				Eigen::Vector3i newFace = Eigen::Vector3i(q_j + V1.rows(), p_i_plus_1,p_i);

				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));
				i = (i + 1) % numBoundaryVerticesScan1;
				cout << "iter = " << iter << endl;
				iter++;
			} while(!INTERP_SURF::edgesAreEqual(newEdge,seedEdge));
		}

		std::cout << "max iter val = " << iter << std::endl;

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

		// assert that indices, used to construct face data, are correct
		for(int i = 0; i < fOff.rows(); i++)
		{
			for(int j = 0; j < 3; j++)
			{
				assert(fOff(i,j) >= 0);
				assert(fOff(i,j) < vOff.rows());
			}	
		}
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
		return (edgeNormDist == 0);
	}


}

﻿// call this <stitching>, not <interp_surface> generation. DEBUG this!
// #TODO :: 
	// [1] actually enable asserts. 
	// [2] Heavily refactor. At least algo still works as expected. I'd refactor it according to pipeline design, but it's kinda too late. Going to have to write <#@$#@> instead atm!
// PURPOSE :: a greedy interpolating surface construction algorithm, from 2 range images, based of sticthing the best edge in a 4-point processing method. #TODO :: write up with a better description!

// MY LIBRARIES
#include "stitchedSurface.h"
#include "remesh.h"
#include "glob_defs.h"

// Other namespace includes
using namespace Eigen;  
using namespace igl;
using namespace std;

// #TODO :: include method to get back boundary vertices
// call is peformed in :: igl::cat(1,bndVertsScan1, bndVertsScan2,vOff);

namespace STITCHED_SURF 
{

	int mod(int a, int b)
	{
		return (a%b+b)%b;
	}

/*
	void generateBoundaryVertices(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, 
								Eigen::MatrixXd& V_boundary)
	{
		Eigen::VectorXi bndIndexes;
		igl::boundary_loop(F,bndIndexes); 
		igl::slice(V,bndIndexes,1,V_boundary); 
	}	
*/
	// uhh, this is a one liner.
	/*
	void retIndicesBoundaryVertices(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, std::vector<int>& bndIndices)
	{
		igl::boundary_loop(F,bndIndices);
	}	
	*/

	void generateStitchedSurface(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1, 
								const Eigen::MatrixXd &V2, const Eigen::MatrixXi &F2, 
								Eigen::MatrixXd& vOff, Eigen::MatrixXi& fOff)
	{
		///////// PREALLOCATION ///////////////////
		struct Mesh
		{
			Eigen::MatrixXd V; 
			Eigen::MatrixXi F;
		} scans;

		std::vector<std::vector<int>> Adjacency_Scan1;
		std::vector<std::vector<int>> Adjacency_Scan2;
		std::vector<bool> boundaryVerticesStatus_scan1;
		std::vector<bool> boundaryVerticesStatus_scan2;

		Eigen::MatrixXd boundaryVertices_scan1;
		Eigen::MatrixXd boundaryVertices_scan2;
		std::vector<int> newTriangleFaces;

	 /*
	  * Interpolating surface generation is a marching face-reconstruction based-approach
	  * we handle bad figure-8 cases, in which one marches error 
	  * this visited_set, is needed only when building a traversing back to the seed edge in 		
	  */

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

		//std::cout << "[1] SOLVE for the seed edge.\n";
	
		std::vector<Eigen::Vector2i> visited;
		Eigen::Vector2i seedEdgeIndices;
		STITCHED_SURF::constructSeedEdgeIndices( bndIndexesScan1, bndVertsScan1, 
							bndIndexesScan2, bndVertsScan2, seedEdgeIndices); 

		Eigen:Vector2i seedEdge;
		//std::cout << "PRE PUSHBACK: " << seedEdge(0) << " " << seedEdge(1) << std::endl;
		//visited.push_back(seedEdge);
		int i = seedEdgeIndices(0); 
		int j = seedEdgeIndices(1);
		int qOffset = numBoundaryVerticesScan1;
		//std::cout << "qOffset = [" << qOffset << "]\n";

		// GOD.FUCKING.DAMMIT.
		// visited.push_back(seedEdge);
		seedEdge(0) = i;
		seedEdge(1) = j + qOffset;
		visited.push_back(seedEdge);
		//std::cout << "POST OFFSET  : " << seedEdge(0) << " " << seedEdge(1) << std::endl;
		Eigen::Vector2i newEdge;
		Eigen::Vector3i newFace; 			

		// std::cout << "[2] Constructed seed Edge Indices" << std::endl;
		// [2] keep alternating edge solving, using the adj lists
		// NOTE :: you DO NOT walk in the same direction, for these cases!
		//	std::cout << "PROGRESSING over Greedy Zippering Surface Reconstruction Algorithm.\n ";

		// std::cout << "E : " << seedEdge(0) << " " << seedEdge(1) << std::endl;
		bool edgesToScan1 = false;
		bool edgesToScan2 = false; 

		// Handle common/seed edge case with two faces
		// ---> contributes to an xWing case, unable to determine orientation
		int pVisit = 0;
		int qVisit = 0;
		do
		{
			// [1] GATHER INDEX AND VERTEX DATA 
			// set up {p_i,p_i+1, q_i,q_i+1}
			int p_i = i;
			int q_j = j;

			int p_i_plus_1  = (STITCHED_SURF::mod(i+1,numBoundaryVerticesScan1));
			int q_j_minus_1 = (STITCHED_SURF::mod(j-1,numBoundaryVerticesScan2));

			//  assess d(b_p_i, b_p_i_plus_1) <= d(b_q_j, b_q_j_minus_1)
			Eigen::VectorXd b_p_i = bndVertsScan1.row(p_i);
			Eigen::VectorXd b_q_j = bndVertsScan2.row(q_j);

			Eigen::VectorXd b_p_i_plus_1 = bndVertsScan1.row(p_i_plus_1);
			Eigen::VectorXd b_q_j_minus_1 = bndVertsScan2.row(q_j_minus_1);
	
			// adding the edges, based on either
			// [a] shortest distance of the diagonal edge
			// [b] whether said edge is a visited edge or not
			double scan1Distance =  (b_q_j - b_p_i_plus_1).squaredNorm(); 
			double scan2Distance =  (b_p_i - b_q_j_minus_1).squaredNorm(); 

			Eigen::Vector2i edgeToScanOne = Eigen::Vector2i(p_i_plus_1, q_j + qOffset);
			Eigen::Vector2i edgeToScanTwo = Eigen::Vector2i(p_i, (q_j_minus_1 + qOffset));

			bool isItEdgeToScanOne  = (scan1Distance < scan2Distance);
			for(std::vector<Eigen::Vector2i>::iterator it = visited.begin(); it != visited.end(); ++it)
			{
				Eigen::Vector2i myEdge = *it;	
				if(edgeToScanOne == myEdge)
				{
					isItEdgeToScanOne = false;
					//std::cout << "Using edge to scan one , repeated case" << std::endl;
					break;
				}
				if(edgeToScanTwo == myEdge)
				{
					isItEdgeToScanOne = true;
					//std::cout << "Using edge to scan two , repeated case" << std::endl;
					break;
				}
			}

			// the bug occurs whenever I flip the triangle faces that I'm working with!
			// relative ordering, of how faces are added, then index check is performed, might matter!
			// #TODO :: why does this NOT execute!!
			if (isItEdgeToScanOne) 
			{
				// newEdge = edgeToScanOne;
			//	std::cout << "Adding edge to Scan 1" << std::endl;
				newEdge = Eigen::Vector2i(p_i_plus_1, q_j + qOffset);
				newFace = Eigen::Vector3i(q_j + qOffset, p_i_plus_1,p_i);

				// check here, if you are back @ seed edge! must be here!
				if (p_i == seedEdge(0)) pVisit++;
				if (pVisit == 2) 
				{
					edgesToScan2 = true;
					//std::cout << "Met same p_i seed edge point" << std::endl;
					break;
				}

				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));
				i = STITCHED_SURF::mod(i + 1, numBoundaryVerticesScan1);

			}
			else 
			{  
				// newEdge = edgeToScanTwo;
				// std::cout << "Adding edge to Scan 2" << std::endl;
				newEdge = Eigen::Vector2i(p_i, (q_j_minus_1 + qOffset));
				newFace = Eigen::Vector3i(p_i,(q_j + qOffset),(q_j_minus_1 + qOffset));

				//std::cout << "q_j = [" << q_j + qOffset << "]\n";
				if (q_j + qOffset == seedEdge(1)) qVisit++;
				if (qVisit == 2) 
				{
					edgesToScan1 = true;
					//std::cout << "Met same q_j seed edge point" << std::endl;
					break;
				}

				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));
				j = STITCHED_SURF::mod(j - 1, numBoundaryVerticesScan2);
			}
			//std::cout << "E : " << newEdge(0) << " " << newEdge(1) << std::endl;
			visited.push_back(newEdge);
		} 
		while(!STITCHED_SURF::edgesAreEqual(newEdge,seedEdge));

		//std::cout << "pVist = [" << pVisit << "]\n.";
		//std::cout << "qVist = [" << qVisit << "]\n.";

		// HACK to dumb corner case ( where you traverse one mesh fully ... b4 the other )
				

		// [2] CLEAN UP PHASE - if the new Edge happens to be adjacent to a seed vertex
		// do not push back 2 repeated triangles --- CGAL will complain!
		if(edgesToScan2)
		{
		//	std::cout << "Cleaning up to Scan 2" << '\n';
			assert(i == seedEdge(0));
			int p_i = i;

			// #TODO :: check if do-while, mustbe converted to a while loop
			// I am basically overwriting my first triangle here ... that is the issue!
			do
			{
				int q_j = j;
				int q_j_minus_1 = (STITCHED_SURF::mod(j-1,numBoundaryVerticesScan2));

				newEdge = Eigen::Vector2i(p_i,(q_j_minus_1 + qOffset));
			//	std::cout << "E : " << newEdge(0) << " " << newEdge(1) << std::endl;
				newFace = Eigen::Vector3i(p_i,(q_j + qOffset),(q_j_minus_1 + qOffset));

				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));
				j = STITCHED_SURF::mod(j - 1, numBoundaryVerticesScan2);
			} 
			while(!STITCHED_SURF::edgesAreEqual(newEdge,seedEdge));
		}
		if (edgesToScan1)
		{
		//	std::cout << "Cleaning up to Scan 1" << '\n';
			assert(bndIndexesScan2(j) == seedEdge(1));
			int q_j = j;

			do
			{
				int p_i = i;
				int p_i_plus_1 = (STITCHED_SURF::mod((i+1),numBoundaryVerticesScan1));

				newEdge = Eigen::Vector2i(p_i_plus_1, (q_j + qOffset));
			//	std::cout << "E : " << newEdge(0) << " " << newEdge(1) << std::endl;
				newFace = Eigen::Vector3i(q_j + qOffset,  p_i_plus_1,p_i);

				// push new face, perform next iteration
				newTriangleFaces.push_back(newFace(0));
				newTriangleFaces.push_back(newFace(1));
				newTriangleFaces.push_back(newFace(2));
				i = STITCHED_SURF::mod(i + 1, numBoundaryVerticesScan1);
			} while(!STITCHED_SURF::edgesAreEqual(newEdge,seedEdge));
		}

		//////////////////////////////////////////////////////////////////////////
		// [3] end once you have the original edge data ! Add this last face  ///
		//////////////////////////////////////////////////////////////////////////
		// CONVERT set of ( 3 * faces ) integers, of vertex indices, to a matrix ( for faces data )

	//	std::cout << "Constructing interpolating surface vertex and face data.\n";
		int numOfFaces = newTriangleFaces.size() / 3;
		Eigen::MatrixXi faces = Eigen::Map<Eigen::MatrixXi,RowMajor> (&newTriangleFaces[0],3,numOfFaces); 

		// CREATE (V,F) off mesh to represent interpolating 
		// surface of the two partial scans.
		igl::cat(1,bndVertsScan1, bndVertsScan2,vOff);
		fOff = faces.transpose();

		// assert that indices, used to construct face data, are correct
		for(int i = 0; i < fOff.rows(); i++)
		{
			for(int j = 0; j < 3; j++)
			{
				assert(fOff(i,j) >= 0);
				assert(fOff(i,j) < vOff.rows());
			}	
		}

		// assert that the orientation is consistent, across the interpolating surface
		for(int i = 0; i < fOff.rows(); ++i)
		{
			for(int j = i + 1; j < fOff.rows(); ++j) 
			{
				Eigen::MatrixXi faceOne = fOff.row(i);
				Eigen::MatrixXi faceTwo = fOff.row(j);

				// check for common subset of two indices
				// a.k.a. you have only sets of size 2, in both cases
	
				// then, based on 2nd common set ... check the ordering! if they're equal ... this is bad!

				std::set<double> setOne;
				std::set<double> setTwo;
				std::set<double>::iterator it;
				std::set<double>::iterator it_i;
				std::set<double>::iterator it_j;

				for(int k = 0; k < 3; ++k)
					setOne.insert(faceOne(k));

				for(int k = 0; k < 3; ++k)
				{
					double curVal = faceTwo(k);
					it = setOne.find(curVal);
					if(it != setOne.end())
						setTwo.insert(curVal);
				}

				if(setTwo.size() == 3)
				{
					std::cout << "REPEATED FACE ERR" << std::endl;
					std::cout << faceOne << std::endl;
					std::cout << faceTwo << std::endl;
				}

				else if (setTwo.size() == 2)
				{
					//std::cout << "COMMON EDGE FOR FACES" << std::endl;
					std::vector<double> commonOne;
					std::vector<double> commonTwo;

					for(int k = 0; k < 3; ++k)
					{
						int v_i = faceOne(k); 
						int v_j = faceOne((k+1) % 3);

						it_i = setTwo.find(v_i);
						it_j = setTwo.find(v_j);
						if(it_i != setTwo.end() && it_j != setTwo.end())
						{
							commonOne.push_back(v_i);
							commonOne.push_back(v_j);
						}
					}

					for(int k = 0; k < 3; ++k)
					{
						int v_i = faceTwo(k); 
						int v_j = faceTwo((k+1) % 3);

						it_i = setTwo.find(v_i);
						it_j = setTwo.find(v_j);
						if(it_i != setTwo.end() && it_j != setTwo.end())
						{
							commonTwo.push_back(v_i);
							commonTwo.push_back(v_j);
						}
					}

					if(commonOne[0] == commonTwo[0] || commonOne[1] == commonTwo[1])
					{
						std::cout << "ERR! Inconsistent orientation amongst faces that share an edge" << std::endl;
						std::cout << faceOne << std::endl;
						std::cout << faceTwo << std::endl;
					}
				}
			}
		}
	
	}

	/*
     * method does not account for offset by vertices in a scan.
     * It assumes that for both scans, the boundary indices 
     * commence from 0 to scan length - 1
     * RESULT :: indices into 'bndVertsSrc' and 'bndVertsDest'
     */
	void constructSeedEdgeIndices(Eigen::VectorXi& bndIndexesSrc, Eigen::MatrixXd& bndVertsSrc,
							Eigen::VectorXi& bndIndexesDest, Eigen::MatrixXd& bndVertsDest,	
							Eigen::Vector2i& seedEdge)
	{
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

/*
 * EXTRANEUOUS NOTES 
 * <igl/combine.h> --- maybe for mesh concatenation. but for now, <igl/cat.h> works.
 * <igl/cut_mesh.h> --- use to cut a mesh? Doesn't really work
 * <igl/orient_outward.h. --- orient mesh normals in a different direction 
 * Are slice methods limited to matrices? Something about tetrahedral meshes too?
 */


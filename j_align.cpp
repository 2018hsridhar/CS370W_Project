// C++ includes
#include <algorithm> 
#include <map>
#include <iostream>

// LibIgl includes
#include <igl/per_face_normals.h>

// MY LIBRARIES  
#include "glob_defs.h"
#include "rigidbodyinstance.h"
#include "j_align.h"
#include "vectormath.h"

// Namespace includes 
using namespace Eigen;  
using namespace std;
using namespace igl;

// PREALLOCATION
//const double h = 0.01;
const double h = 0.1;
//const double h = 0.1;

// NOTE that the <body> is assoc with a <boundary> in this set up.
namespace J_ALIGN
{
	void solve_config_deltas(const Eigen::MatrixXd& vBoundary, const Eigen::MatrixXd bndExtForceMat, const RigidBodyInstance* body, Eigen::Vector3d& delta_cvel, Eigen::Vector3d& delta_omega)
	{
		Eigen::MatrixXd bndExtImpMat = h * bndExtForceMat;
		int numBV = bndExtImpMat.rows();
	
		for(int k = 0; k < numBV; ++k)
		{
			const Eigen::Vector3d J_k_ext = bndExtImpMat.row(k);

			// ADD UP CONTRIBUTIONS TO DELTA_CVEL
			Eigen::Vector3d delta_cVel_k = J_k_ext;
			delta_cvel += delta_cVel_k;
		
			// ADD UP CONTRIBUTIONS TO DELTA_OMEGA
			const Eigen::Vector3d c = vBoundary.row(k); 
			const Eigen::Vector3d C = body->c; 
			const Eigen::Vector3d theta = body->theta;

			Eigen::Vector3d delta_omega_k = VectorMath::rotationMatrix(theta) * ((VectorMath::rotationMatrix(-1.0 * theta) * (C-c)).cross(VectorMath::rotationMatrix(-1.0 * theta) * J_k_ext));
			delta_omega += delta_omega_k;
		}
	}

	// UPDATE CONFIGURATION 
	void updateConfiguration(RigidBodyInstance* ScanBody, const Eigen::Vector3d& delta_cVel,const Eigen::Vector3d& delta_omega)
	{
		ScanBody->cvel += delta_cVel;
		ScanBody->c += (ScanBody->cvel) * h; 
		ScanBody->cvel.setZero();

		// VANISH OUT ( LINEAR & ANGULAR ) VELOCITIES - akin to 'molasses'
		ScanBody->w += delta_omega; 
		const Eigen::Matrix3d rotMat = VectorMath::rotationMatrix(h * ScanBody->w) * VectorMath::rotationMatrix(ScanBody->theta);
		ScanBody->theta = VectorMath::axisAngle(rotMat);
		ScanBody->w.setZero();
	}

	void solveTransformation(const RigidBodyInstance* ScanBody, Eigen::Matrix4d& T)
	{
		const Eigen::Vector3d t_comp = ScanBody->c - ScanBody->c_0;
		const Eigen::Vector3d theta_diff = (ScanBody->theta - ScanBody->theta_0);
		//const Eigen::Vector3d zero_diff = Eigen::Vector3d(0,0,0); // useful for debugging mostly
		//const Eigen::Matrix3d r_comp = VectorMath::rotationMatrix(zero_diff);
		const Eigen::Matrix3d r_comp = VectorMath::rotationMatrix(theta_diff);
		// I expect r_comp = I_3. if not, something else is wrong in the computational of these rotation matrices
		T.block<3,1>(0,3) = t_comp; 
		T.block<3,3>(0,0) = r_comp;
	}

	/*
	*  Just a basic modulus operation
	*/
	int mod(int a, int b)
	{
		return (a%b+b)%b;
	}

	/*
	* GIVEN a specific edge, <e_i> and set of faces <F>
	* RETURN unique face INDEX <f_i> corresponding to those two vertices indexed by the edge <e_i>
	* NOTE that <e_i> is a boundary edge, and thus always indexes into boundary vertices.
	* USE face_index due to <per_face_normals.h> library.
	*/
	int detectFace(const Eigen::Vector2i& e_i, const Eigen::MatrixXi& F)
	{
		// SOLVE for the intersection of a 3d and a 2d vector.
		// Just use a HM freq map . Iterate over <e_i> as if those were your keys instead.
		// COMPLEXITY :: [O(F*E) time, O(E) space]
		bool hasFace = false;
		int f_i = -1; 
		for(int i = 0; i < F.rows(); ++i)
		{
			// CONSTRUCT HM frequency of [face,edge] vertex indices.
			Eigen::Vector3i myFace = F.row(i);
			std::map<int,int> vFreq;
			vFreq[myFace(0)]++;
			vFreq[myFace(1)]++;
			vFreq[myFace(2)]++;
			vFreq[e_i(0)]++;
			vFreq[e_i(1)]++;

			// IF we hit two vertices, we have a matching face
			if(vFreq[e_i(0)] == 2 && vFreq[e_i(1)] == 2) 
			{
				f_i = i;
				break;
			}
		}
		return f_i;
	}

	/*
	* RETURNS a set of boundary vertices and assoc external force normals for each vertex
	* INPUT = [remeshBoundaryOne, remeshed.V, remeshed.F]
	* OUTPUT = [boundryV, extBndForceMat]
	*/
	void solveExternalForcesForBoundary(const std::vector<int>& boundaryLoop, const Eigen::MatrixXd& remeshV, const Eigen::MatrixXi& remeshF, Eigen::MatrixXd& boundaryV, Eigen::MatrixXd& extBndForceMat)
	{
		// [1] Calculate all per face normals, and normalize them
		// #TODO :: <Z> := should I be concerned about cases of faces with degenerate normals?
		// #TODO :: why in the example is it that degenerate faces are assigned (1/3. 1/3, 1/3)^0.5? 
		Eigen::Vector3d Z = Eigen::Vector3d(1,1,1).normalized();
		Eigen::MatrixXd N_face;
		igl::per_face_normals(remeshV, remeshF, Z, N_face);

		// [2] construct set of edges from <boundary_loop_vertices>
		// Eigen::Matrix3i <vFIndices> is of form (v_i, v_j, f_i). It'll serve as a useful index map.
		// #TODO :: ASSERT all faces are unique [ do this later ]
		int numBV = boundaryLoop.size();
		Eigen::MatrixXi vfIndices;	 // size = |bndryLoop| 
		vfIndices.resize(numBV,3);	
		for ( int i = 0; i < numBV; ++i)
		{
			int v_i = boundaryLoop[mod(i,numBV)];
			int v_i_plus_1 = boundaryLoop[mod((i+1), numBV)];
			Eigen::Vector2i e_i = Eigen::Vector2i(v_i,v_i_plus_1);
			int f_i = detectFace(e_i, remeshF);
			if(f_i == -1)
				std::cout << "ERROR :: INVALID FACE INDEX\n";
			Eigen::Vector3i myData = Eigen::Vector3i(v_i,v_i_plus_1,f_i);
			vfIndices.row(i) = myData;
		}

		// [3] Construct normalized edge vectors and face normals.
		// TAKE their cross product :: create new matrix of f_i_ext. 
		// #TODO :: ordering of X-prod vectors is tied to how the boundary is iterated over too.
		Eigen::MatrixXd f_i_ext_Mat;
		Eigen::MatrixXd V_norm;

		f_i_ext_Mat.resize(numBV, 3);
		V_norm.resize(numBV, 3);

		for(int i = 0; i < numBV; ++i)
		{
			Eigen::Vector3i myData = vfIndices.row(i);
			int v_i = myData(0);
			int v_i_plus_1 = myData(1);
			int f_i = myData(2);

			// Solve for face normal
			Eigen::Vector3d f_norm;
			if(f_i == -1) 
			{ 
				// NOTE :: I do not expect this err to arise at all.
				std::cout << "f_i = [-1] ERROR" << std::endl;
				f_norm = Eigen::Vector3d(1,0,0);
				exit(0);
			}
			f_norm = N_face.row(f_i);	

			// solve for edge normal
			Eigen::Vector3d e_i = (remeshV.row(v_i_plus_1) - remeshV.row(v_i));
			double e_i_len = e_i.norm();
			e_i.normalize();

			// Solve for the external force for bndry vertex
			Eigen::Vector3d f_i_ext = (f_norm.cross(e_i)).normalized() * e_i_len; 
			f_i_ext_Mat.row(i) = f_i_ext;

			// Construct boundary vertex assoc w/specific external force <f_i_ext>  
			Eigen::Vector3d c = 0.5 * (remeshV.row(v_i) + remeshV.row(v_i_plus_1));
			V_norm.row(i) =  c;
		}

		// ASSERT :: cardinality(N_edges) matches cardinality(boundary_one_vertices)
		if(numBV != f_i_ext_Mat.rows())
		{
			std::cout << "SIZE MISMATCH ERROR" << std::endl;
			exit(0);
		}

		boundaryV = V_norm;
		extBndForceMat = f_i_ext_Mat;
	}
}

// #TODO :: try to exploit innate properties of structures / utilize acceleration structures to quicken computation of values  ... this takes too damn long to compute !
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include <igl/rotation_matrix_from_directions.h>
#include <igl/writeOFF.h> 
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h> 
#include <igl/procrustes.h>
#include <igl/all_pairs_distances.h> 
#include "helpers.h"
//#include "glob_defs.h"
#include "icp.h"

using namespace Eigen;  
using namespace std;

/*

Eigen::MatrixXd V_one;
	Eigen::MatrixXi F_one;

	Eigen::MatrixXd V_two;
	Eigen::MatrixXi F_two;

	Eigen::MatrixXd V_approx;
	Eigen::MatrixXi F_approx;


bool key_down( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout << "Key : " << key << (unsigned int) key << std::endl;
  if ( key == '1' )
  {
    // clear data before drawing mesh
    viewer.data.clear();
    viewer.data.set_mesh(V_one,F_one);
    viewer.core.align_camera_center(V_one,F_one); 
  }
  else if ( key == '2' ) 
  {
    // clear data before drawing mesh
    viewer.data.clear();
    viewer.data.set_mesh(V_two,F_two);
    viewer.core.align_camera_center(V_two,F_two);
  } 
  else if ( key == '3' ) 
  {
    // clear data before drawing mesh
    viewer.data.clear();
    viewer.data.set_mesh(V_approx,F_approx);
    viewer.core.align_camera_center(V_approx,F_approx);
  }
  return false;
}
*/

/*
int main(int argc, char *argv[])
{

  // PRINT out critical information to user / developer 
  // LOAD mesh data , in OFF format
  std::cout << R"(
1 switch to identity view
2 Switch to rotated view
3 Switch to ICP view 
    )";

  if(!igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V_one, F_one)) {
    cout << " Failed to load mesh [bunny]" << endl;
  }
  if(!igl::readOFF(TUTORIAL_SHARED_PATH "/dataset90Deg.off", V_two, F_two)) {
    cout << " Failed to load mesh [2] " << endl;
  }

  ICP::applyRigidICP(V_one, V_two, V_approx);

  // SETUP LibIgl Viewer 
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data.set_mesh(V_one, F_one);
  viewer.launch();
}
*/

namespace ICP
{
	void applyRigidICP(Eigen::MatrixXd& V_one, Eigen::MatrixXd& V_two, Eigen::MatrixXd& V_transformed)
	{
	
	/***********************************************************/ 
	// PREALLOCATE matrices for computation ( p_i, q, transformation_iter_k)
	/***********************************************************/ 
	
	/* TODO : IS it here that I need to perform barycentric based random sampling ? */
	Eigen::MatrixXd p_i = HELPER::convertToHomogenousForm(V_one); 
	Eigen::MatrixXd q = HELPER::convertToHomogenousForm(V_two); 

	int maxIters = 1000; 										// max num of iteration  
	double tolerance = 0.005;	 								// accuracy threshold  
	Eigen::MatrixXd T_0 = Eigen::MatrixXd::Identity(4, 4);  	// initial guess 
	Eigen::MatrixXd T_k = T_0;								// guess at iter k
	Eigen::VectorXd zeroTranslate = Eigen::Vector4d (0,0,0,0);
	double e_k = 1000; 										// energy to be minimized
	double e_kPlus1 = 1000;								    // energy to be minimized 
	double delta_e = 1000; 								// delta enegy 

	int k = 0;
	Eigen::Matrix4d T = Eigen::MatrixXd::Identity(4,4);

	/***********************************************************/ 
	// ITERATIVELY improve rigid ICP method , for solving the transformation matrix
	/***********************************************************/ 
	while ( k < maxIters && delta_e > tolerance ) 
	{
		std::cout << "Executing iteration " << k << "\n.";
		std::cout << "Delta energy = " << delta_e << "\n.";

		/***********************************************************/ 
		// DISCOVER closest points in q, from p_k 
		/***********************************************************/ 
		Eigen::MatrixXd p_k = (p_i * T_k).rowwise() + zeroTranslate.transpose(); 
		Eigen::MatrixXd p_k_nonHomo = HELPER::normalizeHomogenousMatrix(p_k);
		Eigen::MatrixXd q_nonHomo = HELPER::normalizeHomogenousMatrix(q);
		Eigen::VectorXi Ele = Eigen::VectorXi::LinSpaced(q.rows(),0,q.rows() - 1);
		Eigen::VectorXd smallestSquaredDists;
		Eigen::VectorXi smallestDistIndxs;

		Eigen::MatrixXd q_j_k; 				// this coresponds to closestPointsTo_p_i_From_q ;
		igl::point_mesh_squared_distance(p_k_nonHomo,q_nonHomo ,Ele,smallestSquaredDists,smallestDistIndxs,q_j_k);   // note :: this is T^(k-1)! 
		Eigen::MatrixXd q_j_k_homog = HELPER::convertToHomogenousForm(q_j_k);

		/***********************************************************/ 
		// APPLY Procrstues to solve for T , 
		//     that minimizes norm (T*p_i - q_j_k ) = (p_k - q_j_k)
		// UPDATE transformation matrix via T_k = T * T_{k-1} 
		/***********************************************************/ 

		Eigen::Matrix4d Rotate;
		Eigen::Vector4d Translate;
		double Scale;
		igl::procrustes(p_k, q_j_k_homog, false,false,Scale,Rotate,Translate ); 

		T = Rotate;
		T.col(3) += Translate;  
		Eigen::MatrixXd T_kPlus1 = T * T_k;  

		// CALCULATE energy_kPlus1  AND delta_energy ( for convergence tests ) 
		Eigen::MatrixXd p_kPlus1 = (p_i * T_kPlus1).rowwise() + zeroTranslate.transpose(); 
		Eigen::MatrixXd allPairDistances;
		igl::all_pairs_distances(p_kPlus1, q_j_k_homog, true, allPairDistances); 
		e_kPlus1 = allPairDistances.trace(); 				

		if ( k == 0 ) {
			delta_e = 1000; 
		} else {
			delta_e = std::abs(e_kPlus1 - e_k );
			e_k = e_kPlus1; 
		}     

		T_k = T_kPlus1;
		k = k + 1;
		}

		/***********************************************************/ 
		// CALCULATE the mesh based on RIGID ICP transformation 
		/***********************************************************/ 
	Eigen::MatrixXd rigidIcpMesh = (p_i * T_k).rowwise() + zeroTranslate.transpose(); 
	V_transformed = HELPER::normalizeHomogenousMatrix(rigidIcpMesh);
	}
}

/***********************************************************/ 
// #TODO :: APPROACH 2 ( using surface normals ... code this up later ) 
//Eigen::MatrixXd meshOneVertexNormals_prime = transformMat_k * meshOneVertexNormals;
// find intersection qi_k on surface Q, based on normal lines !
// STEP (1) :: SELECT control points pi \in P ( i = 1..N), , compute surface normals n_pi, set initial transformation matrix T_0
// I am very unsure? Plus there is the issue of getting normals too , from barycentric ( not sure about that too ! ) 
/***********************************************************/ 


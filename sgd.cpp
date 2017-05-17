// Q1 :: What is the order of convergence of this method? not sure atm
// purpose of this code :: given a (4x4) transformation matrix/alignemnt, 
// iteratively improve the alignemnt, via stoichastic gradient descent !  
// try 100 of these things ... take min. it will most likely converge ! 

#include "glob_defs.h"
#include "sgd.h"
#include <Eigen/Geometry>

using namespace Eigen; 
using namespace std;

//const int numTransMats = 5;
const int numTransMats = 15;

namespace SGD
{
	void generateTransMats(Eigen::Matrix4d& input, std::vector<Eigen::Matrix4d>& transMats)
	{
		// SET UP RNG , for epsilon values {x,y,z,\theta}
		// - [-0.001, 0.001] too low
		// - [-1, 1] too high 
		// - [-0.1,0.1] ... maybe okay? some explosion?
		// - [-0.01, 0.01] ... to low?
		//const double range_from  = -0.03;
		//const double range_to    =  0.03;
		const double range_from  = -0.05;
		const double range_to    =  0.05;
		std::random_device rand_dev;
		std::mt19937 generator(rand_dev());
		std::uniform_real_distribution<double>  distr(range_from, range_to);

		for(int k = 0; k < numTransMats; k++)
		{
			// GENERATE 3 values for epsilon distances
			// - three values for {x,y,z} (rotation and translation)
			// - one value for an epsilon angle ( rotation ) 
			double e_x = distr(generator);
			double e_y = distr(generator);
			double e_z = distr(generator);
			//double epsilon_axis[] = { e_x,e_y,e_z};
			//Eigen::Vector3d e_vec(e_x,e_y,e_z);
			double e_theta = distr(generator); 

			// CREATE rotation component, 
			// [axis-angle] -> [quaternion] -> [homogenous transformation]

			// check for NUMERICAL ISSUES !! ... where are you getting suspicious numbers!
			// - check rotation matrices are orthogonal!
			// - to use valgrind on this?? ... hmmm ... 
 
			//Eigen::Vector3d epsilon_axis2 = igl::random_dir(); 
/*
			double qrot[4];
			igl::axis_angle_to_quat(epsilon_axis,epsilon_angle,qrot); // does libIGl expect a positive angle? 

			// print out qrot
			cout << "qrot coeff vals= " << endl;
			for(int i = 0; i < 4; ++i)
				cout << qrot[i] << ',';
			cout << endl;
			double mat[16];
			igl::quat_to_mat(qrot, mat);
			cout << "Rot coeff mats = " << endl;
			for(int i = 0; i < 16; ++i)
				cout << mat[i] << ',';
			cout << endl;
*/

			Eigen::Matrix4d rotation = Eigen::Matrix4d::Identity();  
			// #TODO :: use segment here!
/*
			for ( unsigned i = 0; i < 4; ++i )
			{
				for ( unsigned j = 0; j < 4; ++j )
				{
					rotation(i,j) = mat[i+(4*j)]; 
				}
			}
*/

			// try a diff way of generation this ... use Eigen's AngleAxis formulae
			Eigen::Matrix3d m;
			Eigen::Vector3d sA(e_x,e_y,e_z);
			sA.normalize();
			m = Eigen::AngleAxisd(e_theta * 10 * igl::PI, sA).toRotationMatrix();

			for ( unsigned i = 0; i < 3; ++i )
			{
				for ( unsigned j = 0; j < 3; ++j )
				{
					rotation(i,j) = m(i,j);
				}
			}

	

			// MAKE translation component
			Eigen::Matrix4d translation = Eigen::Matrix4d::Zero();
			translation(0,3) = e_x;
			translation(1,3) = e_y;
			translation(2,3) = e_z;

			// ASSERT that rotation and translations make sense
/*
			Eigen::Matrix4d assertOrtho = rotation.transpose() * rotation;
			std::cout << "rotation matrix = " << std::endl;
			std::cout << rotation << std::endl;
			std::cout << "Orthog rotation = " << std::endl;
			std::cout << assertOrtho << std::endl;
			std::cout << "translation matrix = " << std::endl;
			std::cout << translation << std::endl;
*/
			// GENERATE perturbed matrix 
			// - perturbed by rotation and translation
			Eigen::Matrix4d perturbed =  (input * rotation) + translation;
 			Eigen::Matrix4d perturbedByR = input * rotation;
        	Eigen::Matrix4d perturbedByT = input + translation;
//			cout << perturbedByT << endl;

			transMats.push_back(perturbed);
			transMats.push_back(perturbedByR);
			transMats.push_back(perturbedByT);
		}
	}

	// ENERGY is calculated merely off of surface area 
	// - can extend further, if need be 
	// - possibly useful libraries ( for later ) ?
	// 		- principal_curvature.h
	// 		- per_corner/edge/face/vertex normals
	// 		- doublearea.h
	//		- gaussian_curvature.h
	// 		- note :: if you don't return anything ... c++ chooses to return u default val ( here, 0 ) LOL
	double calculateSurfaceEnergy(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
	{
		// REFERENCE THEIR KNOWN BUG ( scales V + F, not just F )
		// solve for surface area of input mesh
		Eigen::VectorXd dbla;
		igl::doublearea(V,F,dbla);
		double surfaceArea = 0.5 * dbla.sum();	

		// calculat energy 
		double result = surfaceArea;
		assert(result >= 0);
		return result;
	}

	// given 2 vectors ( transformation matrices + energies), 
	// discover optimal transformation matrix ( and return it's energy val, as a double )

    // BE REALLY CAREFUL CODING THIS UP!!!
	double findOptimalTransMat(std::vector<Eigen::Matrix4d>& transMats, std::vector<double>& energies, Eigen::Matrix4d& opt)
	{
		assert(transMats.size() == energies.size());

		// print out the transMats and energy values
/*
		for(int i = 0; i < energies.size(); ++i)
		{
			std::cout << "i = [" << i << "]" << std::endl;
			std::cout << transMats[i] << std::endl;
			std::cout << energies[i] << std::endl;
			std::cout << "---------------------------------" << std::endl;
		}
*/
		// get index of max energy value
		//int idxMinEnergy = -1; 
		int idxMinEnergy = 0; 
		double minEnergy = std::numeric_limits<double>::max();
		for(int i = 0; i < energies.size(); ++i)
		{
			if(energies[i] < minEnergy)
			{
				minEnergy = energies[i];
				idxMinEnergy = i;
			}
		}

		assert(idxMinEnergy >= 0);
		assert(idxMinEnergy < energies.size());

		opt = transMats[idxMinEnergy];
		return minEnergy;
	}
}

// REFERNCE FOR RANDOMIZER :: 
// http://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range/20136256#20136256
// generate N random (3x1) vector of small {e_x,e_y,e_z} values
// note to self ... iterators suck for vectors WHEN you need to index into the vectors themselves! 

// good tip ... if a function returns some notion of an error ( i.e. could not execute a call ... capture that). JUST notion of not being able to execute the call though!

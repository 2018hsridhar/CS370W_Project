// Q1 :: What is the order of convergence of this method? not sure atm
// purpose of this code :: given a (4x4) transformation matrix/alignemnt, 
// iteratively improve the alignemnt, via stoichastic gradient descent !  
// try 100 of these things ... take min. it will most likely converge ! 

#include "glob_defs.h"
#include "sgd.h"

using namespace Eigen; 
using namespace std;

const int numTransMats = 10;
//const int numTransMats = 4;
namespace SGD
{
	void generateTransMats(Eigen::Matrix4d& input, std::vector<Eigen::Matrix4d>& transMats)
	{
		// SET UP RNG , for epsilon values {x,y,z,\theta}
		// - [-0.001, 0.001] too low
		// - [-1, 1] too high 
		const double range_from  = -0.1;
		const double range_to    =  0.1;
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
			double epsilon_axis[] = { e_x,e_y,e_z};
			double epsilon_angle = distr(generator); 

			// CREATE rotation component, 
			// [axis-angle] -> [quaternion] -> [homogenous transformation]

			// check for NUMERICAL ISSUES !! ... where are you getting suspicious numbers!
			// - check rotation matrices are orthogonal!
			// - to use valgrind on this?? ... hmmm ... 
 
			Eigen::Vector3d epsilon_axis2 = igl::random_dir(); 
			double qrot[4];
			igl::axis_angle_to_quat(epsilon_axis,epsilon_angle,qrot);
			double mat[16];
			igl::quat_to_mat(qrot, mat);

			Eigen::Matrix4d rotation = Eigen::Matrix4d::Identity();  
			for ( unsigned i = 0; i < 4; ++i )
			{
				for ( unsigned j = 0; j < 4; ++j )
				{
					rotation(i,j) = mat[i+(4*j)]; 
				}
			}

			// MAKE translation component
			Eigen::Matrix4d translation = Eigen::Matrix4d::Zero();
			translation(0,3) = e_x;
			translation(1,3) = e_y;
			translation(2,3) = e_z;

			// ASSERT that rotation and translations make sense
			Eigen::Matrix4d assertOrtho = rotation.transpose() * rotation;
			std::cout << "Orthog rotation = " << std::endl;
			std::cout << assertOrtho << std::endl;
			std::cout << "translation matrix = " << std::endl;
			std::cout << translation << std::endl;
		
			// GENERATE perturbed matrix 
			// - perturbed by rotation and translation
			Eigen::Matrix4d perturbed =  (input * rotation) + translation;
			transMats.push_back(perturbed);
		}
	}

	// ENERGY is calculated merely off of surface area 
	// - can extend further, if need be 
	// - possibly useful libraries ( for later ) ?
	// 		- principal_curvature.h
	// 		- per_corner/edge/face/vertex normals
	// 		- doublearea.h
	//		- gaussian_curvature.h
	double calculateSurfaceEnergy(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
	{
		// REFERENCE THEIR KNOWN BUG ( scales V + F, not just F )
		// solve for surface area of input mesh
		Eigen::VectorXd dbla;
		igl::doublearea(V,F,dbla);
		// note :: does this return an error code !
		double surfaceArea = 0.5 * dbla.sum();	

		// calculat energy 
		double result = -1;
		result = surfaceArea;
		assert(result >= 0);
		return result;
	}

	// given 2 vectors ( transformation matrices + energies), discover optimal transformation matrix
	void findOptimalTransMat(std::vector<Eigen::Matrix4d>& transMats, std::vector<double>& energies, Eigen::Matrix4d& opt)
	{
		assert(transMats.size() == energies.size());

		// get index of max energy value
		int idxMinEnergy = -1;
		double minEnergy = std::numeric_limits<double>::max();
		for(int i = 0; i < energies.size(); ++i)
		{
			if(energies[i] < minEnergy)
			{
				minEnergy = energies[i];
				idxMinEnergy = i;
			}
		}
		assert(idxMinEnergy >= 0 && idxMinEnergy < energies.size());

		// return corresponding permutation matrix
		opt = transMats[idxMinEnergy];	
	}
}

// REFERNCE FOR RANDOMIZER :: 
// http://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range/20136256#20136256
// generate N random (3x1) vector of small {e_x,e_y,e_z} values
// note to self ... iterators suck for vectors WHEN you need to index into the vectors themselves! 

// good tip ... if a function returns some notion of an error ( i.e. could not execute a call ... capture that). JUST notion of not being able to execute the call though!

// Q1 :: What is the order of convergence of this method? not sure atm
// purpose of this code :: given a (4x4) transformation matrix/alignemnt, 
// iteratively improve the alignemnt, via stoichastic gradient descent !  
// try 100 of these things ... take min. it will most likely converge ! 

#include "glob_defs.h"
#include "sgd.h"

using namespace Eigen; 
using namespace std;

//const int numTransMats = 10;
const int numTransMats = 4;
namespace SGD
{
	void generateTransMats(Eigen::Matrix4d& input, std::vector<Eigen::Matrix4d>& transMats)
	{
		// SET UP random number generator , for epsilon values
		//const double range_from  = -0.001;
		//const double range_to    = 0.001;
		const double range_from  = -0.1;
		const double range_to    = 0.1;
		//const double range_from  = -1; // TO DAMM HIGIH!
		//const double range_to    = 1;
		std::random_device rand_dev;
		std::mt19937 generator(rand_dev());
		std::uniform_real_distribution<double>  distr(range_from, range_to);

		// #TODO :: refactor code, to be cleaner ... ehhh, do that later
// of course you will have an isuse here, since glob_defs.h does not get linked to this program, ONLY to the main program!
		//for(int i = 0; i < GLOBAL::numTransMats; i++)
		for(int i = 0; i < numTransMats; i++)
		{
			// generate 3 values for epsilon distances
			// in {x,y,z} directions
			// and an epsilon angle too
			double e_x = distr(generator);
			double e_y = distr(generator);
			double e_z = distr(generator);
			double epsilon_axis[] = { e_x,e_y,e_z};
			double epsilon_angle = distr(generator); 

			// generate rotation component, via axis-angle and quat conversions
			Eigen::Vector3d epsilon_axis2 = igl::random_dir(); 
// this might be a better idea ! and then convert to an arr?? 
			double qrot[4];
			igl::axis_angle_to_quat(epsilon_axis,epsilon_angle,qrot);
			double mat[16];
			igl::quat_to_mat(qrot, mat);

			Eigen::Matrix4d rotatePerturbation = Eigen::Matrix4d::Identity();  // note :: translate vector is zeroed out here! 
			for ( unsigned i = 0; i < 4; ++i )
				for ( unsigned j = 0; j < 4; ++j )
				rotatePerturbation(i,j) = mat[i+4*j]; 

			// generate translation permutation component
			Eigen::Matrix4d translatePerturbation = Eigen::Matrix4d::Zero();
			translatePerturbation(0,3) = e_x;
			translatePerturbation(1,3) = e_y;
			translatePerturbation(2,3) = e_z;

			// 
		
			// generate perturbed matrix ( perturbed by rotation and translation )
			Eigen::Matrix4d perturbedRotAndTransMat =  (input * rotatePerturbation) + translatePerturbation;
			Eigen::Matrix4d assertOrtho = rotatePerturbation.transpose() * rotatePerturbation;
			std::cout << "Orthog rotation = " << std::endl;
			std::cout << assertOrtho << std::endl;
			transMats.push_back(perturbedRotAndTransMat);
		}
	}

	// ENERGY is calculated merely off of surface area - can extend further, if need be #TODO
	// possibly useful libraries ( for later ) ?
	// 		- principal_curvature.h
	// 		- per_corner/edge/face/vertex normals
	// 		- doublearea.h
	//		- gaussian_curvature.h
	double calculateSurfaceEnergy(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
	{
		// solve for surface area of input mesh
		Eigen::VectorXd dbla;
		igl::doublearea(V,F,dbla);
		double surfaceArea = 0.5 * dbla.sum();	

		// calculat energy 
		double result = -1;
		result = surfaceArea;
		assert(result >= 0);
		return result;
	}

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

// REFERNCE THIS FOR RANDOMIZER :: 
// http://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range/20136256#20136256
// generate N random (3x1) vector of small {e_x,e_y,e_z} values
// note to self ... iterators suck for vectors WHEN you need to index into the vectors themselves! 


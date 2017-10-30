#include "helpers.h"
#include <Eigen/Core>
#include <Eigen/Geometry> 

namespace HELPER 
{
	Eigen::MatrixXd convertToHomogenousForm(const Ref<const MatrixXd>& mat)
	{
		Eigen::MatrixXd homogenoizedMatrix;
		homogenoizedMatrix = mat.transpose().colwise().homogeneous().transpose(); 
		return homogenoizedMatrix;
	}

	Eigen::MatrixXd normalizeHomogenousMatrix (const Ref<const MatrixXd>& mat)
	{
		Eigen::MatrixXd normalizedMatrix;
		normalizedMatrix = mat.transpose().colwise().hnormalized().transpose(); 
		return normalizedMatrix;
	}

	void applyRigidTransformation(const Eigen::MatrixXd& V, const Eigen::MatrixXd& T, Eigen::MatrixXd& vOut)
	{
		Eigen::VectorXd zeroTranslate = Eigen::Vector4d (0,0,0,0);
		Eigen::MatrixXd p_i = HELPER::convertToHomogenousForm(V);
		//std::cout << "converted to homogenous form is" << std::endl;
		//std::cout << p_i << std::endl;
		//Eigen::MatrixXd transV = (T * p_i).rowwise() + zeroTranslate.transpose(); // THIS IS WRONG!
		Eigen::MatrixXd transV = (T * p_i.transpose()).transpose(); // THIS IS RIGHT!!
		//Eigen::MatrixXd transV = p_i;
		vOut = HELPER::normalizeHomogenousMatrix(transV);
	//	std::cout << "converted to regular form is " << std::endl;
	//	std::cout << vOut << std::endl;
	}

	void constructNormalField(Eigen::MatrixXd& V, Eigen::MatrixXd& N, Eigen::MatrixXd& V_res, Eigen::MatrixXi& N_edges)
	{
		// #TODO :: convert this to it's own method body


		// #NOTE :: you'll convert this into a method later that takes in (V,N) and returns (V',E) to represent these visualized normals.
		// apply the normal scaling operation here
		Eigen::MatrixXd offsetedV = V + N;
		//scan2.V = (scan1.V).rowwise() + n.transpose(); 

		// concat the two vertex meshes too ( need this for indexing to properly work )
		igl::cat(1,V,offsetedV,V_res);

		// set up edges between [Scan1.V,scan2.V]
		//std::cout << "COMMENCING set up of edges corresponding to normals\n";
		int offset = V.rows();
		N_edges = Eigen::MatrixXi(offset,3);
		for(int i = 0; i < offset; ++i)
		{
			Eigen::VectorXi e_i = Eigen::Vector3i(i,i + offset, i);  // #TODO :: fix this problem
			// hmm ... cheat, set one of these the same ????
			N_edges.row(i) = e_i;	
		}	
		//std::cout << "DONE setting up edges corresponding to normals\n";
	}

}

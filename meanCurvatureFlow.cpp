// STANDARD INCLUDE FILES 
#include "meanCurvatureFlow.h"
#include "tutorial_shared_path.h"
#include <igl/is_border_vertex.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h> 

using namespace Eigen; 
using namespace std;

namespace MCF
{
	bool meshHasBoundary(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
	{
		std::vector<bool> boundary = igl::is_border_vertex(V, F);  
		// NEED counter here, as boundary isa  boolean vector ACROSS all vertices
		int boundCount = 0;
		for(int i = 0; i < boundary.size(); i++)
		{
			if(boundary[i])
			{
				boundCount++;
			}
		}
		return ( boundCount > 0);
	}

	void computeMeanCurvatureFlowWaterTight(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double timestep, Eigen::MatrixXd &Vc)
	{
		// PREALLOCATION
		Eigen::SparseMatrix<double> M, L;
		igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
		Vc = V;
		double diff;
		constexpr double mcfStoppingVal = 32;  
		// #TODO :: try to choose a better convergence/stopping value, over here!

		//CALCUALTE the initial mass and stiffness matrices
		igl::cotmatrix(V,F,L);  
		igl::massmatrix(V,F, mcfType, M);  

		// ITERATE till reach convergence criteria
		do
		{
			diff = 0;

			// SOLVE for new vertices, U, based on linear, discretized equation 
			Eigen::SparseMatrix<double> A = ( M - ( timestep * L));
			Eigen::MatrixXd B = ( M * Vc); 
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver; 
			solver.compute(A); 
			if(solver.info() != Eigen::Success) {
				std::cout << "Decomposition of A failed." << std::endl;
			}
			auto updatedMeshVertices = solver.solve(B);
			if ( solver.info() != Eigen::Success )  {
				std::cout << "Solving B failed." << std::endl;
			}
			auto U = updatedMeshVertices.eval(); 
			if ( solver.info() != Eigen::Success )  {
				std::cout << "Solving B failed." << std::endl;
			}

			// RESCALE U, by previous unit-area
			Eigen::VectorXd dbla;
			igl::doublearea(U,F,dbla);
			double area = 0.5 * dbla.sum(); 
			U /= std::sqrt(area);

			// SOLVE for difference, between current and previous iteration of MCF ( Vc vs U )
			// NOTE :: does difference make sense ... get to a really tiny, singular val here!
			int numVertices = Vc.rows();
			for(int i = 0; i < numVertices; i++)
			{
				diff += (Vc.row(i) - U.row(i)).norm();		
			}
			printf("Difference this round is %f.\n", diff);
			Vc = U;

			// RECALCULATE mass matrix
			igl::massmatrix(Vc,F, mcfType, M);  
		} while ( diff > mcfStoppingVal);

	}

	void computeMeanCurvatureFlowBoundary(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, std::vector<bool> boundary, double timestep, Eigen::MatrixXd &Vc)
	{
		// PREALLOCATION
		Eigen::SparseMatrix<double> M, L;
		igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
		Vc = V;
		double diff;
		constexpr double mcfStoppingVal = 10;  
		// #TODO :: fix difference of norm calculation ... it becomes constant, @ some point of time ... this is not good !

		//CALCUALTE the initial mass and stiffness matrices
		igl::cotmatrix(V,F,L);  
		igl::massmatrix(V,F, mcfType, M);  

		// ITERATE till reach convergence criteria
		do
		{
			diff = 0;
			Eigen::SparseMatrix<double> A = ( M - ( timestep * L));
			Eigen::MatrixXd B = ( M * Vc); 
			Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver; 
			// #TODO :: why does Clements use a SparseLU solver, specifically? 
			solver.compute(A); 
			if(solver.info() != Eigen::Success) {
				std::cout << "Decomposition of A failed." << std::endl;
			}
			// #TODO :: DISCOVER difference in these 
			// 		two if-conditinos ( solve and evaluate stage ) ! 
			auto updatedMeshVertices = solver.solve(B);
			if ( solver.info() != Eigen::Success )  {
				std::cout << "Solving B failed." << std::endl;
			}
			auto U = updatedMeshVertices.eval(); 
			if ( solver.info() != Eigen::Success )  {
				std::cout << "Solving B failed." << std::endl;
			}
			int numVertices = Vc.rows();
			// assert(solver.info() == Eigen::Success) 
			// #TODO :: figure out how to use "assert()"

			// SOLVE for difference, between current and previous iteration of MCF ( Vc vs U )
			for(int i = 0; i < numVertices; i++)
			{
				diff += (Vc.row(i) - U.row(i)).norm();		
			}
			printf("Difference this round is %f.\n", diff);
			
			// UPDATE Vc = (U | Vc ), based on if a vertex is on the boundary
			for(int i = 0; i < numVertices; i++) 
			{
				if ( !boundary[i] ) 
				{
					Vc.row(i) = U.row(i);
				}
			}

			// RECALCULATE mass matrix
			igl::massmatrix(Vc,F, mcfType, M);  
		} while ( diff > mcfStoppingVal);
	}

	void computeMeanCurvatureFlow(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F, double timestep, Eigen::MatrixXd &Vc)
	{
		bool hasBndryVertices = meshHasBoundary(V,F);
		if(hasBndryVertices) {
		    std::vector<bool> boundary = igl::is_border_vertex(V,F);	
			computeMeanCurvatureFlowBoundary(V,F,boundary,timestep,Vc);
		}
		else { 
			computeMeanCurvatureFlowWaterTight(V,F,timestep,Vc);
		}
	}
}
	/* 
	****************************** 
	***** BELONGS ELSEWHERE *****
	****************************** 
	*/
	// STL caauses issues ... OBJ and OFF formats do not ... not sure why though !
	//igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off", V_one, F_one); 
	//igl::readOFF(TUTORIAL_SHARED_PATH "/proper_sphere.off", V_one, F_one);  // it straight up does not work for this case !
	//igl::readSTL("decimated-max.stl", V_one, F_one,N_one);  // it straight up does not work for this case !
	//igl::readOBJ(TUTORIAL_SHARED_PATH "/decimated-max.obj",V_one,F_one);
	//igl::readOFF("horseAhead.off", V_one, F_one);  // it straight up does not work for this case !
	//igl::readOFF(TUTORIAL_SHARED_PATH "/beetle.off",V_one,F_one); /// this finally works 


/*
	... OLD VOL RESCALING TECHNIQUE ... 

			// update vertices ( coefficent vector ) and mass matrix
			Eigen::VectorXd dbla;
			igl::doublearea(V_mcf,F_one,dbla);
			double oldVolume = 0.5 * dbla.sum(); 

			V_mcf = newVertices;
			igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
			igl::massmatrix(V_mcf,F_one, mcfType, massMatrix_iterK);  
			//std::cout << "[3] Succesfully updated mass and stiffness matrices \n";

			Eigen::VectorXd dbla_new;
			igl::doublearea(V_mcf,F_one,dbla_new);
			double newVolume = 0.5 * dbla_new.sum(); 

			V_mcf /= std::cbrt(newVolume); 
*/

// new change in design :: the user should be able to apply mean curvature flow, but not neccessarily based on whether there is or there is not a boundary! hence, hide those as private methods in the namespace!  ... #TODO :: get to thsi change ,,, eventually !


/*
	bool key_down( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
	{
	std::cout << "Key : " << key << (unsigned int) key << std::endl;
	if ( key == '1' ) {
		viewer.data.clear();
		viewer.data.set_mesh(V_one,F_one);
		viewer.core.align_camera_center(V_one,F_one); 
	}
	else if ( key == '2' ) {
		viewer.data.clear();
		//applyOneTimeStepOfMcfWatertightCase();
		applyOneTimeStepOfMcfBoundaryCase();
		viewer.data.set_mesh(V_mcf,F_one);
		viewer.core.align_camera_center(V_mcf,F_one);
		std::string output_of_mesh = "meanCurvaureFlowOutput[" + std::to_string(k) + "].stl";
		std::cout << "WRITING MCF data for iteration [" << std::to_string(k) << " ].\n"; 
		igl::writeOFF(output_of_mesh, V_mcf, F_one);
	} 
	else if ( key == '3' ) {
		viewer.data.clear();
		V_mcf = V_one;
		k = 0;
		igl::cotmatrix(V_one,F_one, stiffnessMatrix_iterK );  
		igl::MassMatrixType mcfType = igl::MASSMATRIX_TYPE_BARYCENTRIC;
		igl::massmatrix(V_one,F_one, mcfType, massMatrix_iterK);  
		viewer.data.set_mesh(V_one,F_one);
		viewer.core.align_camera_center(V_mcf,F_one);
	} 
	return false;
	}
	*/
	
	/*
	* REFACTORING CODE HERE ! 
	*
	*/




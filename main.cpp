﻿// MY LIBRARIES  
#include "glob_defs.h"
#include "meanCurvatureFlow.h"
#include "interpSurface.h"
#include "remesh.h"
#include "sgd.h"
#include "helpers.h"

// LibIgl includes
#include <igl/writeOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/readOFF.h>

using namespace Eigen;  
using namespace std;
using namespace igl;

struct Mesh
{
	Eigen::MatrixXd V; 
	Eigen::MatrixXi F;
} scan1,scan2,interp,remeshed, result;

// testing one iteration of SGD - from 1 matrices ( init @ I ) to 5 matrices ( @ what SGD produes for thee ). I should desend specifically on the one with LOWEST energy!
void runPipeline();
//bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier);

int main(int argc, char *argv[])
{
	runPipeline();
	
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;	
	igl::readOFF(GLOBAL::pipelineOutputFile, V, F);
	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(V,F);
	viewer.launch();
}

/*

// in a couple of wyas, setting up this visualizer seems like a dumb idea
// it is low-hanging fruit, but  ... IT IS ALSO NOT THAT CRUCIAL!!

bool key_down( igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	std::cout << "Key : " << key << (unsigned int) key << std::endl;
	if ( key == 'r' || key == 'R' ) 
	{
		viewer.data.clear();
	}
	else if ( key == ' ' )
	{
		viewer.data.clear();
		applyOneTimeStepOfPipeline();
		viewer.data.set_mesh(V_mcf,F_one);
		viewer.core.align_camera_center(V_mcf,F_one);
		std::string output_of_mesh = "meanCurvaureFlowOutput[" + std::to_string(k) + "].stl";
		std::cout << "WRITING MCF data for iteration [" << std::to_string(k) << " ].\n"; 
		igl::writeOFF(output_of_mesh, V_mcf, F_one);
	} 
	return false;
}
*/
	
// let me modify the method, to view stepping - BUT, don't PRINT out the output file, just yet!
void runPipeline()
{
	std::cout << "Executing Zero-Overlap Rigid Alignment Pipeline" << std::endl;

	if(!readOFF(GLOBAL::pipelineScan1File,scan1.V,scan1.F)) {
		std::cout << "Failed to load partial scan one." << std::endl;
	} 

	if(!readOFF(GLOBAL::pipelineScan2File,scan2.V,scan2.F)) {
		std::cout << "Failed to load partial scan two." << std::endl;
	}

	std::cout << "Executing pipeline for the following meshes" << endl;
	std::cout << "Mesh one [" << GLOBAL::pipelineScan1File << "]" << std::endl;
	std::cout << "Mesh two [" << GLOBAL::pipelineScan2File << "]" << std::endl;

	// Develop an initial transformation matrix, filled with random data 
	Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
	std::vector<Eigen::Matrix4d> transMats;	
	transMats.push_back(T);

	// k = 10, 15 ... get reasonable output
	// k = 20 = GAAH
	for(int k = 0; k < 20; ++k)
	{
		cout << "Running SGD for [" << k << "] th iteration" << endl;;
		std::vector<double> energies;
		for(int i = 0; i < transMats.size(); ++i)
		{
			cout << "i = [" << i << "]" << endl;
			cout << "Applying pipeline, with transition matrix" << endl;
			cout << transMats[i] << endl;

			cout << "Apply Rigid Transformation to Scan 1 Mesh" << endl;
			Eigen::MatrixXd transScan1;	
			HELPER::applyRigidTransformation(scan1.V,T,transScan1);

			cout << "Generating interpolating surface" << endl;
			INTERP_SURF::generateOffsetSurface(transScan1,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
			cout << "Remeshing interpolating surface" << endl;
			double remeshEdgeLen = REMESH::avgEdgeLenInputMeshes(transScan1,scan1.F,scan2.V,scan2.F);
			REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, remeshEdgeLen);

			cout << "Flowing the interpolating surface" << endl;
			Eigen::MatrixXd V = remeshed.V;
			Eigen::MatrixXi F = remeshed.F;
			bool hasBndry = MCF::meshHasBoundary(V,F);
			assert(hasBndry);
			Eigen::MatrixXd Vc;
			MCF::computeMeanCurvatureFlow(V,F,0.001, Vc);
		
			cout << "Calculating energy value" << endl;
			double energy = SGD::calculateSurfaceEnergy(Vc,F);
			energies.push_back(energy);
		}
		SGD::findOptimalTransMat(transMats,energies,T);
		transMats.clear();
		SGD::generateTransMats(T,transMats);
		cout << "OPTIMAL transition matrix is : " << endl;
		cout << T << endl;
		cout << "-----------------------------------------------" << endl;
	}
	cout << "Pipeline execution finished" << endl;	
	cout << "END RESULT transition matrix is : " << endl;
	cout << T << endl;

	// write data to output file  - create one huge mesh
	// result.V, result.F
	Eigen::MatrixXd transScan1;	
	HELPER::applyRigidTransformation(scan1.V,T,transScan1);
	igl::cat(1,transScan1,scan2.V,result.V);
	igl::cat(1,scan1.F, MatrixXi(scan2.F.array() + scan1.V.rows()), result.F);
	igl::writeOFF(GLOBAL::pipelineOutputFile, result.V, result.F);

	cout << "Wrote off file for aligned meshes" << endl;
}



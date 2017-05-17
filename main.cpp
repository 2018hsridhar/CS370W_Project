﻿// reference for c++ timer library is :: http://www.cplusplus.com/reference/ctime/
// http://www.cplusplus.com/reference/ctime/clock/ ... uhh, make sure your units make sense! 
// C++ includes
#include <iostream>
#include <fstream>
#include <time.h>

// MY LIBRARIES  
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

void runPipeline();

int main(int argc, char *argv[])
{
	runPipeline();
}

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
	double prev_energy = std::numeric_limits<double>::max();
	int k = 0;

	// change animated pipeline to break, based on whether local energy is a minimum!
	// REMESHING seems to be SUCH bottleneck here!
	clock_t init, final;
	init=clock();
	while(true)
	{
		std::vector<double> energies;
		for(int i = 0; i < transMats.size(); ++i)
		{
			Eigen::MatrixXd transScan1;	
			HELPER::applyRigidTransformation(scan1.V,transMats[i],transScan1);
			INTERP_SURF::generateOffsetSurface(transScan1,scan1.F,scan2.V,scan2.F, interp.V,interp.F);
			double remeshEdgeLen = REMESH::avgEdgeLenInputMeshes(transScan1,scan1.F,scan2.V,scan2.F);
			REMESH::remeshSurface(interp.V,interp.F,remeshed.V,remeshed.F, remeshEdgeLen);

			// let's just assume we already have a boundary in these cases ! 
			Eigen::MatrixXd Vc;
			MCF::computeMeanCurvatureFlow(remeshed.V,remeshed.F,0.001, Vc);
		
			double energy = SGD::calculateSurfaceEnergy(Vc,remeshed.F);
			energies.push_back(energy);
		}
		Eigen::Matrix4d T_local;
		double local_energy = SGD::findOptimalTransMat(transMats,energies,T_local);
		if(local_energy < prev_energy)
		{
			T = T_local;
			prev_energy = local_energy;
		}
		else
		{
			// you found your local matrix here! 
			break;
		}
		transMats.clear();
		SGD::generateTransMats(T,transMats);
		energies.clear();
		++k;
	}
	final=clock()-init;
	double totalExecTime = (double)final / ((double)CLOCKS_PER_SEC);
	
	std::cout << "Pipeline execution finished.\n";
	std::cout << "END RESULT transition matrix is : \n";
	std::cout << T << std::endl;
	printf("Total Number Iterations = [%d]\n",k);
	printf("Total Execution time = [%f] seconds.\n", totalExecTime);
}

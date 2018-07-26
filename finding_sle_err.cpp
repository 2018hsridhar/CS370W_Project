﻿// CODE PURPOSES :: emulate CGAL's SLE operation, for simple triangle mesh
// --- CGAL introduces possible FLE issues here
// --- Can vary edge length accordingly.
// EXEC on top of one simple test case focus only [ 1 triangle case :: fine-tune half-scale]
// INDEX faces properly
// EXECUTE a testing barrage 
// good tests files = [2triangles.off, fandisk.off, planexy.off]. 
#include "glob_defs.h"
#include "helpers.h"
#include "stitchedSurface.h"
#include "remesh.h"
#include "meanCurvatureFlow.h"
#include "sgd.h"  
#include "vectormath.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "j_align.h"

// viewer logic ( viewing mostly )
#include <igl/viewer/Viewer.h>
#include <igl/cat.h>
#include <igl/boundary_loop.h>

// LibIgl includes [ specific for this code too ! ] 
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/avg_edge_length.h>

// LibIgl - nice printing statements too!
#include <igl/print_ijv.h>
#include <igl/print_vector.h> // wait a sec, this library looks incorrect [ refer to commentary ] 

// C++ includes
#include <iostream>
#include <fstream>
#include <algorithm> 

// C++ Namespace includes 
using namespace Eigen;  
using namespace std;
using namespace igl;

// Forward declaration of classes

// PREALLOCATION
struct Mesh
{
	Eigen::MatrixXd V; 
	Eigen::MatrixXi F;
	Eigen::MatrixXd bV;  // boundary vertices
	Eigen::VectorXi bI; // boundary indices 
	Eigen::MatrixXd N;
	int numBV; 			// number of boundary vertices
};
struct Mesh myImg; 

// METHOD HEADERS
int runPipeline();
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int mod);
ofstream debugFile;
int mod(int a, int b);

// SIMPLE MOD OPERATION
int mod(int a, int b)
{
	return (a%b+b)%b;
}

// SET 2 options - [1] for default view, [2] for first SLE view, [3] for the second SLE view.
// #TODO :: find out more about <OFF> file metadata. I do not think this matters as much as I thought.
namespace SLE
{
	void execSLE()
	{
		// LOAD the simple triangle mesh in and view it.
		if(!readOFF(GLOBAL::simpleTriangleFile,myImg.V, myImg.F)) {
			std::cout << "Failed to load myImg." << std::endl;
			return -1;
		}
		readOFF(GLOBAL::simpleTriangleFile,myImg.V,myImg.F);
		std::cout << "Executing basic SLE tests for <myImg>" << std::endl;

		// SOLVE for the average edge length over all edges
		double rEL = igl::avg_edge_length(myImg.V, myImg.F);
		double k = 2.0; // magnitutde to divide <rEL> by here
		rEL = rEL / k; // scale down edge length by 1/2

		// ADD vertices in <SLE> operations
		// NOTE --- you must create a new set of {V,F} here now, and we do not know the exact size of {V,F} here too. This might present some future issues! Uhh, isn't it just the scalar I divide by from <avg_edge_length> calculations? ... hmmm ...?
		// better idea - just create 2 dynamically updated array lists for <V,F> respectively. That would prove far easier to process here!
		// ... more ever, for SLE operation, as you create a vertex, you really (1) delete one triangle and (2) create 2 new triangles in place. But in the next step, you (1) delete that triagnel, and create (2) new ones in place again. This algorithm is rather inductive, isn't it?
		// MORE-EVER, I can be given an arbitrary mesh [ e.g. planexy ], so I need to focus specifically on boundary-vertices only!  ... but how do we get that other vertex [ the interior vertex, specifically ]. I have to iterate over boundary_faces somehow too. #TODO circumnavig this hurdle too.
		// another issue ... how to properly iterate in the triangles case? I expect [1->6] faces for a [1/2] SLE split functions, right?? 
		// NOTE - we can safely ignore issues regarding the proper aspect ratio. This really is not as problematic as we expected. We still have remeshing available
		// NOTE that edge lengths VARY - so the #-vertices, per edge, will also accordingly vary
		// YES, the behavior from the [single_triangle] case is as expected - let's write that up! 


		// <---------------------------------------------------------------------------->
		// NOTE that we first need to identify the <boundary_loops> for a given [V,E] mesh. Iterate ACCORDING to that!!
		// ... BUT how to properly iterate over <boundary_edges>? That is an issue. BUT ... <boundary_loop>'s themselves provide indexes themselves into each boundary edge. So that's useful!
		std::vector<int> inBoundryLoop;
		igl::boundary_loop(myImg.F, inBoundryLoop);
	/*
		std::cout << "[";
		for ( int i : inBoundryLoop)
			std::cout << i << ", ";
		std::cout << "]\n";
	*/

		// NOW - loop over edges. Get this edge len. 
		// USE <avg_edge_len> to introduce your vertices.
		// <EDGE_LEN> = (v2-v1).norm()

		// well, we always know to start from v1, so add v1, then the rest, in your list!. 
		std::vector<Eigen::Vector3d> myVerts;
		int num_BV = inBoundryLoop.size();
		for ( int i = 0; i < num_BV; ++i)
		{
			Eigen::Vector3d v1 = myImg.V.row(i);
			Eigen::Vector3d v2 = myImg.V.row(mod(i + 1,num_BV));
			double e_len = (v2-v1).norm();
			int num_verts = e_len / rEL;
			std::cout << num_verts << std::endl; 
			// use n-midpoint rules here, for solving these vertices.
			for(int a = 0; a < num_verts; ++a)
			{
				Eigen::Vector3d offset = ( v2 - v1 ) * ( 1.0 / e_len );
				Eigen::Vector3d v_prime = v1 + a * (offset);
				myVerts.push_back(v_prime);
			}
		} 

		// PRINT OUT current vertex list
		int new_V_card = myVerts.size();
		for(int i = 0; i < new_V_card; ++i)
		{
			Eigen::Vector3d myV = myVerts[i];
			std::cout << "[" << i << "]\t[" << myV(0) << ", " << myV(1) << ", " << myV(2) << "]" << std::endl;
			// CREATE the new face with [vertex - 1, curVertex, vertex + 1]
			// YOU have to do this while in the processing stage. 
		} 

		// PERFORM edge-splitting
		// let us make edge of these edges ... HUH, perhaps NOT!
		// split W.R.T. faces too! 
		// 2 possible techniques atm?
		// 		[1] split faces as vertices are added? 
		// 		[2] add all vertices, THEN split the faces?

		// LAUNCH igl viewer	
		igl::viewer::Viewer viewer;
		viewer.data.clear();
		viewer.data.set_mesh(myImg.V,myImg.F);
		viewer.launch();
	}
}

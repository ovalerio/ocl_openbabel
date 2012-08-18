/*********************************************************************
forcefieldmmff94gpu.cpp - MMFF94 force field using GPU

Based on forcefieldmmff94.cpp
Copyright (C) 2006-2008 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

Modifications to run using GPU acceleration
Copyright (C) 2012 by Omar Valerio <omar.valerio@gmail.com>


This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/*
 * Source code layout:
 * - Functions to calculate the actual interactions
 * - Parse parameter files
 * - Setup Functions
 * - Validation functions
 * - Calculate bond type, angle type, stretch-bend type, torsion type
 * - Various tables & misc. functions
 *
 */

#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/locale.h>

#include <iomanip>
#include <cmath>

#include "forcefieldmmff94gpu.h"

//#define __NO_STD_VECTOR // Use cl::vector instead of STL version
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>

//
// Include Eigen headers
//
#include <Eigen/Dense>

// Benchmark utils contais the class used for timing kernel runs
//#include "benchmark-utils.hpp"

using namespace std;
using namespace cl;

namespace OpenBabel
{
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Functions to calculate the actual interactions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////


  double OBForceFieldMMFF94GPU::Energy(bool gradients)
  {
    double energy = 0;

    IF_OBFF_LOGLVL_MEDIUM {
		snprintf(_logbuf, BUFF_SIZE, "\n USE GRADIENTS = %s \n ",(gradients)?"true":"false");
		OBFFLog(_logbuf);
    }

    IF_OBFF_LOGLVL_MEDIUM {
      OBFFLog("\nE N E R G Y\n\n");
    }

//	energy += EnergyBond();
//	energy += EnergyAngle();  // !!! THIS METHOD IMPLEMENTATION IS NOT READY
//	energy += EnergyStrBnd();
//	energy += EnergyTorsion();
//	energy += EnergyOOP();
//	energy += EnergyVDW();
	energy += EnergyElectrostatic();

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "\nTOTAL ENERGY = %8.5f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //
  // MMFF part I - page 494
  //
  //                   kb_ij                              7
  // EB_ij = 143.9325 ------- /\r_ij^2 (1 + cs /\_rij + ---- cs^2 r_ij^2)
  //                     2                               12
  //
  // kb_ij	force constant (md/A)
  //
  // /\r_ij 	r_ij - r0_ij (A)
  //
  // cs		cubic stretch constant = -2 A^(-1)
  //

  double OBForceFieldMMFF94GPU::EnergyBond(){
	    double energy = 0.0;

	    std::vector<OBFFBondCalculationMMFF94>::iterator iterator;

	    std::vector<Eigen::Vector3d> posAVector;
	    std::vector<Eigen::Vector3d> posBVector;
	    std::vector<double> r0Vector;
	    std::vector<double> kbVector;

	    Eigen::Vector3d posA, posB;

	    for (iterator = _bondcalculations.begin(); iterator != _bondcalculations.end(); iterator++){
	    	posA << (*iterator).pos_a[0], (*iterator).pos_a[1], (*iterator).pos_a[2] ;
	    	posAVector.push_back(posA);
	    	posB << (*iterator).pos_b[0], (*iterator).pos_b[1], (*iterator).pos_b[2] ;
	    	posBVector.push_back(posB);
	    	r0Vector.push_back((*iterator).r0);
	    	kbVector.push_back((*iterator).kb);
	    }//for(bondCalculations)

	    Eigen::Vector3d rAB;

	    double r0, kb;
	    double delta, delta_sq;

	    for (int indx = 0; indx < _bondcalculations.size(); indx++){
	    	posA = posAVector[indx];
	    	posB = posBVector[indx];
	    	r0 = r0Vector[indx];
	    	kb = kbVector[indx];

	    	// calculate the bond vectors and determine its length (distance)
	    	rAB = posA - posB;
	    	delta = rAB.norm() - r0;
	    	delta_sq = pow(delta, 2);

	    	energy += kb * delta_sq * (1.0 - 2.0 * delta + 7.0/3.0 * delta_sq);

	    }//for(bondCalculations)

	    energy = 143.9325 * 0.5 * energy;

	    return energy;

  }

  //
  // MMFF part I - page 495
  //
  //                       ka_ijk
  // EA_ijk = 0.438449325 -------- /\0_ijk^2 (1 + cs /\0_ijk)
  //                         2
  //
  // ka_ijk	force constant (md A/rad^2)
  //
  // /\0_ijk 	0_ijk - 00_ijk (degrees)
  //
  // cs		cubic bend constant = -0.007 deg^-1 = -0.4 rad^-1
  //

  double OBForceFieldMMFF94GPU::EnergyAngle(){

	    double energy = 0.0;

	    std::vector<OBFFAngleCalculationMMFF94>::iterator iterator;

	    std::vector<Eigen::Vector3d> posAVector;
	    std::vector<Eigen::Vector3d> posBVector;
	    std::vector<Eigen::Vector3d> posCVector;
	    std::vector<double> theta0Vector;
	    std::vector<double> kaVector;
	    std::vector<bool> linearVector;

	    Eigen::Vector3d posA, posB, posC;

	    for (iterator = _anglecalculations.begin(); iterator != _anglecalculations.end(); iterator++){
	    	posA << (*iterator).pos_a[0], (*iterator).pos_a[1], (*iterator).pos_a[2] ;
	    	posAVector.push_back(posA);
	    	posB << (*iterator).pos_b[0], (*iterator).pos_b[1], (*iterator).pos_b[2] ;
	    	posBVector.push_back(posB);
	    	posC << (*iterator).pos_c[0], (*iterator).pos_c[1], (*iterator).pos_c[2] ;
	    	posCVector.push_back(posC);
	    	theta0Vector.push_back((*iterator).theta0);
	    	kaVector.push_back((*iterator).ka);
	    	linearVector.push_back((*iterator).linear);
	    }//for(angleCalculations)

	    Eigen::Vector3d rAB, rBC;

	    double theta, cos_theta;
	    double theta0, delta_theta;
	    double ka;
	    bool linear;

	    for (int indx = 0; indx < _anglecalculations.size(); indx++){
	    	posA = posAVector[indx];
	    	posB = posBVector[indx];
	    	posC = posCVector[indx];
	    	theta0 = theta0Vector[indx];
	    	ka = kaVector[indx];
	    	linear = linearVector[indx];

	    	// calculate the bond vectors between the three atoms
	    	rAB = posA - posB;
	    	rBC = posC - posB;

	    	// perform in place normalization of the vectors
	    	rAB.normalize();
	    	rBC.normalize();

	    	// calculate the cosine_theta and obtain theta
	    	cos_theta = rAB.dot(rBC);
	    	theta = RAD_TO_DEG * acos(cos_theta);

	    	delta_theta = theta - theta0;

	        if (linear) {
	        	energy += 143.9325 * ka * (1.0 + cos(theta * DEG_TO_RAD));
	        } else {
		    	energy += 0.043844 * 0.5 * ka * pow(delta_theta, 2) * (1.0 - 0.007 * delta_theta);
	        }

//	        std::cout << "ka:" << ka << "theta0:" << theta0 << "theta:" << theta <<  "linear:" << linear std::endl;

	    }//for(angleCalculations)

	    return energy;

  }

  //
  // MMFF part I - page 495
  //
  // EBA_ijk = 2.51210 (kba_ijk /\r_ij + kba_kji /\r_kj) /\0_ijk
  //
  // kba_ijk	force constant (md/rad)
  // kba_kji	force constant (md/rad)
  //
  // /\r_xx 	see above
  // /\0_ijk 	see above
  //

  double OBForceFieldMMFF94GPU::EnergyStrBnd(){

	    double energy = 0.0;

	    std::vector<OBFFStrBndCalculationMMFF94>::iterator iterator;

	    std::vector<Eigen::Vector3d> posAVector;
	    std::vector<Eigen::Vector3d> posBVector;
	    std::vector<Eigen::Vector3d> posCVector;
	    std::vector<double> theta0Vector;
	    std::vector<double> rab0Vector;
	    std::vector<double> rbc0Vector;
	    std::vector<double> kbaABCVector;
	    std::vector<double> kbaCBAVector;

	    Eigen::Vector3d posA, posB, posC;

	    for (iterator = _strbndcalculations.begin(); iterator != _strbndcalculations.end(); iterator++){
	    	posA << (*iterator).pos_a[0], (*iterator).pos_a[1], (*iterator).pos_a[2] ;
	    	posAVector.push_back(posA);
	    	posB << (*iterator).pos_b[0], (*iterator).pos_b[1], (*iterator).pos_b[2] ;
	    	posBVector.push_back(posB);
	    	posC << (*iterator).pos_c[0], (*iterator).pos_c[1], (*iterator).pos_c[2] ;
	    	posCVector.push_back(posC);
	    	theta0Vector.push_back((*iterator).theta0);
	    	rab0Vector.push_back((*iterator).rab0);
	    	rbc0Vector.push_back((*iterator).rbc0);
	    	kbaABCVector.push_back((*iterator).kbaABC);
	    	kbaCBAVector.push_back((*iterator).kbaCBA);
	    }//for(strbndCalculations)

	    Eigen::Vector3d rAB, rBC;
	    double lenAB, lenBC;

	    double cos_theta, theta;
	    double delta_theta, delta_rab, delta_rbc;
	    double theta0, rab0, rbc0;
	    double factor, kbaABC, kbaCBA;

	    for (int indx = 0; indx < _strbndcalculations.size(); indx++){
	    	posA = posAVector[indx];
	    	posB = posBVector[indx];
	    	posC = posCVector[indx];
	    	theta0 = theta0Vector[indx];
	    	rab0 = rab0Vector[indx];
	    	rbc0 = rbc0Vector[indx];
	    	kbaABC = kbaABCVector[indx];
	    	kbaCBA = kbaCBAVector[indx];

	    	// calculate the bond vectors between the three atoms
	    	rAB = posA - posB;
	    	rBC = posC - posB;

	    	// calculate vector lengths for the bond vectors
	    	lenAB = rAB.norm();
	    	lenBC = rBC.norm();

	    	// perform in place normalization of the vectors
	    	rAB.normalize();
	    	rBC.normalize();

	    	// calculate the cosine_theta and obtain theta
	    	cos_theta = rAB.dot(rBC);
    		theta = RAD_TO_DEG * acos(cos_theta);

	    	delta_theta = theta - theta0;
	    	delta_rab = lenAB - rab0;
	    	delta_rbc = lenBC - rbc0;
	    	factor = kbaABC * delta_rab + kbaCBA * delta_rbc;

//	    	std::cout << "theta:" << theta << "theta0:" << theta0 << "dtheta:" << delta_theta << " dRAB:" << delta_rab << " dRBC" << delta_rbc << std::endl;

	    	energy += factor * delta_theta;

	    }//for(strbndCalculations)

	    energy = 2.51210 * energy;

	    return energy;

  }


  //
  // MMFF part I - page 495
  //
  // ET_ijkl = 0.5 ( V1 (1 + cos(0_ijkl)) + V2 (1 - cos(2 0_ijkl)) + V3 (1 + cos(3 0_ijkl)) )
  //
  // V1		force constant (md/rad)
  // V2		force constant (md/rad)
  // V3		force constant (md/rad)
  //
  // 0_ijkl 	torsion angle (degrees)
  //

  double OBForceFieldMMFF94GPU::EnergyTorsion(){
	    double energy = 0.0;

	    std::vector<OBFFTorsionCalculationMMFF94>::iterator iterator;

	    std::vector<Eigen::Vector3d> posAVector;
	    std::vector<Eigen::Vector3d> posBVector;
	    std::vector<Eigen::Vector3d> posCVector;
	    std::vector<Eigen::Vector3d> posDVector;
	    std::vector<Eigen::Vector3d> velVector;

	    Eigen::Vector3d posA, posB, posC, posD, vel;

	    for (iterator = _torsioncalculations.begin(); iterator != _torsioncalculations.end(); iterator++){
	    	posA << (*iterator).pos_a[0], (*iterator).pos_a[1], (*iterator).pos_a[2] ;
	    	posAVector.push_back(posA);
	    	posB << (*iterator).pos_b[0], (*iterator).pos_b[1], (*iterator).pos_b[2] ;
	    	posBVector.push_back(posB);
	    	posC << (*iterator).pos_c[0], (*iterator).pos_c[1], (*iterator).pos_c[2] ;
	    	posCVector.push_back(posC);
	    	posD << (*iterator).pos_d[0], (*iterator).pos_d[1], (*iterator).pos_d[2] ;
	    	posDVector.push_back(posD);
	    	vel << (*iterator).v1, (*iterator).v2, (*iterator).v3 ;
	    	velVector.push_back(vel);
	    }//for(torsionCalculations)

	    Eigen::Vector3d rAB, rBC, rCD;
	    Eigen::Vector3d nR, nS, nT;
	    double lenAB, lenBC, lenCD;

	    double d1, d2, tor;
	    double cos1, cos2, cos3;
	    Eigen::Vector3d phi;

	    for (int indx = 0; indx < _torsioncalculations.size(); indx++){
	    	posA = posAVector[indx];
	    	posB = posBVector[indx];
	    	posC = posCVector[indx];
	    	posD = posDVector[indx];
	    	vel = velVector[indx];

	    	// calculate bond vectors between the three atoms
	    	rAB = posB - posA;
	    	rBC = posC - posB;
	    	rCD = posD - posC;

	    	// calculate vector lengths for the bond vectors
	    	lenAB = rAB.norm();
	    	lenBC = rBC.norm();
	    	lenCD = rCD.norm();

	    	// perform in place normalization of the vectors
	    	rAB.normalize();
	    	rBC.normalize();
	    	rCD.normalize();

	    	// calculate the normal vectors of the three planes
	    	nR = rAB.cross(rBC);
	    	nS = rBC.cross(rCD);
	    	nT = nR.cross(nS);

	    	// calculate d1, d2, tor
	    	d1 = nT.dot(rBC);
	    	d2 = nR.dot(nS);
	    	tor = RAD_TO_DEG * atan2(d1, d2);

	    	cos1 = cos(DEG_TO_RAD * 1 * tor);
	    	cos2 = cos(DEG_TO_RAD * 2 * tor);
	    	cos3 = cos(DEG_TO_RAD * 3 * tor);

	    	phi(0) = 1.0 + cos1;
	    	phi(1) = 1.0 - cos2;
	    	phi(2) = 1.0 + cos3;

	    	energy += vel.dot(phi);

	    }//for(torsionCalculations)

	    energy = 0.5 * energy;

	    return energy;
  }


  //						//
  //  a						//
  //   \  					//
  //    b---d      plane = a-b-c		//
  //   / 					//
  //  c						//
  //						//

  double OBForceFieldMMFF94GPU::EnergyOOP(){
	  double energy = 0.0;

	    std::vector<OBFFOOPCalculationMMFF94>::iterator iterator;

	    std::vector<Eigen::Vector3d> posAVector;
	    std::vector<Eigen::Vector3d> posBVector;
	    std::vector<Eigen::Vector3d> posCVector;
	    std::vector<Eigen::Vector3d> posDVector;
	    std::vector<double> koopVector;

	    Eigen::Vector3d posA, posB, posC, posD;

	    for (iterator = _oopcalculations.begin(); iterator != _oopcalculations.end(); iterator++){
	    	posA << (*iterator).pos_a[0], (*iterator).pos_a[1], (*iterator).pos_a[2] ;
	    	posAVector.push_back(posA);
	    	posB << (*iterator).pos_b[0], (*iterator).pos_b[1], (*iterator).pos_b[2] ;
	    	posBVector.push_back(posB);
	    	posC << (*iterator).pos_c[0], (*iterator).pos_c[1], (*iterator).pos_c[2] ;
	    	posCVector.push_back(posC);
	    	posD << (*iterator).pos_d[0], (*iterator).pos_d[1], (*iterator).pos_d[2] ;
	    	posDVector.push_back(posD);
	    	koopVector.push_back((*iterator).koop);
	    }

	    Eigen::Vector3d rBA, rBC, rBD;
	    Eigen::Vector3d nABC, nCBD, nABD;
	    double angle, normBA, normBC, normBD;

	    double theta, cos_theta, sin_theta, sin_dl;
	    double koop;

	    // The Out of Plane term compute was originally adapted from Andreas Moll dissertation on BALLView
	    // http://scidok.sulb.uni-saarland.de/volltexte/2007/1325/pdf/Dissertation_1544_Moll_Andr_2007.pdf
	    for (int indx = 0; indx < _oopcalculations.size(); indx++){
	    	posA = posAVector[indx];
	    	posB = posBVector[indx];
	    	posC = posCVector[indx];
	    	posD = posDVector[indx];
	    	koop = koopVector[indx];

	    	// calculate bond vectors from central atom to outer atoms
	    	rBA = posA - posB;
	    	rBC = posC - posB;
	    	rBD = posD - posB;

	    	// calculate vector lengths for the bond vectors
	    	normBA = rBA.norm();
	    	normBC = rBC.norm();
	    	normBD = rBD.norm();

	    	// perform in place normalization of the vectors
	    	rBA.normalize();
	    	rBC.normalize();
	    	rBD.normalize();

	    	// calculate the normal vectors of the three planes
	    	nABC = rBA.cross(rBC);
	    	nCBD = rBC.cross(rBD);
	    	nABD = rBD.cross(rBA);

	    	// theta is the angle between rBA and rBC
	    	cos_theta = rBA.dot(rBC);
	    	theta = acos(cos_theta);
	    	sin_theta = sin(theta);

	    	sin_dl = nABC.dot(rBD) / sin_theta;

	        // the wilson angle is the asin
	        angle = RAD_TO_DEG * asin(sin_dl);

//	      if (!isfinite(angle))
//	        angle = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

	      energy += koop * pow(angle, 2);

	    }// for(vdwcalculations)

	    energy = 0.043844 * 0.5 * energy;

	    return energy;
  }


  /// Van der Waals non-bonded interaction contribution

  double OBForceFieldMMFF94GPU::EnergyVDW(){
	  double energy = 0.0;

	    std::vector<OBFFVDWCalculationMMFF94>::iterator iterator;

	    std::vector<Eigen::Vector3d> posAVector;
	    std::vector<Eigen::Vector3d> posBVector;
	    std::vector<double> RABVector, RAB7Vector;
	    std::vector<double> epsilonVector;

	    Eigen::Vector3d posA, posB;

	    for (iterator = _vdwcalculations.begin(); iterator != _vdwcalculations.end(); iterator++){
	    	posA << (*iterator).pos_a[0], (*iterator).pos_a[1], (*iterator).pos_a[2] ;
	    	posAVector.push_back(posA);
	    	posB << (*iterator).pos_b[0], (*iterator).pos_b[1], (*iterator).pos_b[2] ;
	    	posBVector.push_back(posB);
	    	RABVector.push_back((*iterator).R_AB);
	    	RAB7Vector.push_back((*iterator).R_AB7);
	    	epsilonVector.push_back((*iterator).epsilon);
	    }

	    Eigen::Vector3d diffAB;
	    double RAB, RAB7;
	    double rab, rab7, qq;
	    double erep, erep7, eattr;
	    double epsilon;

	    for (int indx = 0; indx < _vdwcalculations.size(); indx++){
	    	posA = posAVector[indx];
	    	posB = posBVector[indx];
	    	diffAB = posA - posB;
	    	rab = diffAB.norm(); // distance between atoms
	    	rab7 = pow (rab, 7);
	    	RAB = RABVector[indx];
	    	RAB7 = RAB7Vector[indx];
	    	epsilon = epsilonVector[indx];
	    	erep = (1.07 * RAB) / (rab + 0.07 * RAB);
	    	erep7 = pow (erep, 7);
	    	eattr = (((1.12 * RAB7) / (rab7 + 0.12 * RAB7)) - 2.0);
	    	energy += epsilon * erep7 * eattr;
	    }// for(vdwcalculations)

	    return energy;
  }


  // Electrostatic contribution

  double OBForceFieldMMFF94GPU::EnergyElectrostatic(){
	    ScalarType energy = 0.0;

	    std::vector<OBFFElectrostaticCalculationMMFF94>::iterator iterator;

	    const int SIZE = _electrostaticcalculations.size();
	    ScalarType * posA = new ScalarType[SIZE*3];
	    ScalarType * posB = new ScalarType[SIZE*3];
	    ScalarType * qq = new ScalarType[SIZE];

	    int andx, bndx, qndx;
	    for (andx=0, bndx = 0, qndx = 0, iterator = _electrostaticcalculations.begin();
	    		iterator != _electrostaticcalculations.end(); iterator++){
	    	posA[andx++] = (*iterator).pos_a[0];
	    	posA[andx++] = (*iterator).pos_a[1];
	    	posA[andx++] = (*iterator).pos_a[2];
	    	posB[bndx++] = (*iterator).pos_b[0];
	    	posB[bndx++] = (*iterator).pos_b[1];
	    	posB[bndx++] = (*iterator).pos_b[2];
	    	qq[qndx++] = (*iterator).qq;
	    }//for(electrostaticCalculations)

	   try {
	        // Get available platforms
	        std::vector<Platform> platforms;
	        Platform::get(&platforms);

	        // Select the default platform and create a context using this platform and the GPU
	        cl_context_properties cps[3] = {
	            CL_CONTEXT_PLATFORM,
	            (cl_context_properties)(platforms[0])(),
	            0
	        };
	        Context context( CL_DEVICE_TYPE_GPU, cps);

	        // Get a list of devices on this platform
	        std::vector<Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

	        // Create a command queue and use the first device
	        CommandQueue queue = CommandQueue(context, devices[0]);

	        // Read source file
			const char * relative_path = "electrostatic_forces_kernel.cl";
		   char full_path[500]; // probably long enough to hold a full path
		   full_path[0] = '\0'; //strcat searches for '\0' to cat after
//		   strcat( full_path, CL_SOURCE_DIR ); // Copy base path into full path
		   strcat( full_path, "kernels");
//		   strcat( full_path, "/media/radar/git/openbabelgpu/openbabel-2.3.1/src/forcefields");
		   strcat( full_path, "/"); // Assume a backslash is used to separate path names
		   strcat( full_path, relative_path); // Copy the relative path to the full_path
//		   std::cout << "path: " << full_path << std::endl;

	        std::ifstream sourceFile(full_path);
	        std::string sourceCode(
	            std::istreambuf_iterator<char>(sourceFile),
	            (std::istreambuf_iterator<char>()));
	        Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

	        // Make program of the source code in the context
	        Program program = Program(context, source);

	        // Build program for these specific devices
	        program.build(devices);

	        // Make kernel
	        Kernel kernel(program, "forces");

	        // Create memory buffers
	        Buffer bufferPosA = Buffer(context, CL_MEM_READ_ONLY, 3 * SIZE * sizeof(ScalarType));
	        Buffer bufferPosB = Buffer(context, CL_MEM_READ_ONLY, 3 * SIZE * sizeof(ScalarType));
	        Buffer bufferQQ = Buffer(context, CL_MEM_READ_ONLY, SIZE * sizeof(ScalarType));
	        Buffer bufferPartial = Buffer(context, CL_MEM_READ_WRITE, SIZE * sizeof(ScalarType));

	        // Copy posA, posB and QQ to the memory buffers
	        queue.enqueueWriteBuffer(bufferPosA, CL_TRUE, 0, 3 * SIZE * sizeof(ScalarType), posA);
	        queue.enqueueWriteBuffer(bufferPosB, CL_TRUE, 0, 3 * SIZE * sizeof(ScalarType), posB);
	        queue.enqueueWriteBuffer(bufferQQ, CL_TRUE, 0, SIZE * sizeof(ScalarType), qq);

	        // Set arguments to kernel
	        kernel.setArg(0, bufferPosA);
	        kernel.setArg(1, bufferPosB);
	        kernel.setArg(2, bufferQQ);
	        kernel.setArg(3, bufferPartial);

	        // Run the kernel on specific ND range
	        NDRange global(SIZE);
	        NDRange local(1);

	        queue.enqueueNDRangeKernel(kernel, NullRange, global, local);

	        // Read partial energies into a local array
	        ScalarType *partial = new ScalarType[SIZE];
	        queue.enqueueReadBuffer(bufferPartial, CL_TRUE, 0, SIZE * sizeof(ScalarType), partial);

	        for(int i = 0; i < SIZE; i ++){
	        	std::cout << partial[i] << std::endl;
	             energy += partial[i];
	        }

	    } catch(Error error) {
	       std::cout << error.what() << "(" << error.err() << ")" << std::endl;
	       std::cout << "Error: " << oclErrorString(error.err()) << std::endl;
	    }


//	    Eigen::MatrixXd ABMatrix;
//	    Eigen::VectorXd rAB;
//	    Eigen::VectorXd energyVector;
//	    Eigen::VectorXd pointZeroFiveVector  = Eigen::VectorXd::Constant(size, 0.05);
//
//	    // calculate the bond vectors and determine its length (distance)
//	    ABMatrix = posAMatrix - posBMatrix;
//	    rAB = ABMatrix.rowwise().norm();
//	    energyVector = qqVector.array() / (rAB + pointZeroFiveVector).array(); // add 0.05 to avoid zero division
//	    energy = energyVector.sum();

	    return energy;
  }


  //
  // OBForceFieldMMFF member functions
  //
  //***********************************************
  //Make a global instance
  OBForceFieldMMFF94GPU theForceFieldMMFF94GPU("MMFF94GPU", false);
  OBForceFieldMMFF94GPU theForceFieldMMFF94sGPU("MMFF94sGPU", false);
  //***********************************************

  OBForceFieldMMFF94GPU::~OBForceFieldMMFF94GPU()
  {
  }

  OBForceFieldMMFF94GPU &OBForceFieldMMFF94GPU::operator=(OBForceFieldMMFF94GPU &src)
  {
    _mol = src._mol;
    _init = src._init;
    return *this;
  }


} // end namespace OpenBabel

//! \file forcefieldmmff94gpu.cpp
//! \brief MMFF94 force field

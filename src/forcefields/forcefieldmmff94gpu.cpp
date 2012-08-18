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

	energy += EnergyBond();
	energy += EnergyAngle();  // !!! THIS METHOD IMPLEMENTATION IS NOT READY
	energy += EnergyStrBnd();
	energy += EnergyTorsion();
	energy += EnergyOOP();
	energy += EnergyVDW();
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

	      energy += koop * pow(angle, 2);

	    }// for(oopcalculations)

	    energy = 0.043844 * 0.5 * energy;

	    return energy;
  }


  /// Van der Waals non-bonded interaction contribution

  double OBForceFieldMMFF94GPU::EnergyVDW(){
	  double energy = 0.0;

		const int size = vdwCalculations.totalPairs();
		FloatType * ax = &vdwCalculations.ax[0];
		FloatType * ay = &vdwCalculations.ay[0];
		FloatType * az = &vdwCalculations.az[0];
		FloatType * bx = &vdwCalculations.bx[0];
		FloatType * by = &vdwCalculations.by[0];
		FloatType * bz = &vdwCalculations.bz[0];
		FloatType * RAB = &vdwCalculations.RAB[0];
		FloatType * epsilon = &vdwCalculations.epsilon[0];

		FloatType rab, rab7, RAB7;
		FloatType abx, aby, abz;
	    FloatType erep, erep7, eattr;

		for (int indx = 0; indx < size; indx++) {
			abx = ax[indx] - bx[indx];
			aby = ay[indx] - by[indx];
			abz = az[indx] - bz[indx];
			rab = sqrt(abx * abx + aby * aby + abz * abz);
	    	rab7 = pow (rab, 7);
	    	RAB7 = pow(RAB[indx], 7);
	    	erep = (1.07 * RAB[indx]) / (rab + 0.07 * RAB[indx]);
	    	erep7 = pow (erep, 7);
	    	eattr = (((1.12 * RAB7) / (rab7 + 0.12 * RAB7)) - 2.0);
	    	energy += epsilon[indx] * erep7 * eattr;
	    }// for(vdwcalculations)

	    return energy;
  }


  // Electrostatic contribution
  // NOTE: For this part we use conversions from vector to array . For example: int* pv = &v[0];
  // This pointer is only valid as long as the vector is not reallocated. Reallocation happens
  // automatically if you insert more elements than will fit in the vector's remaining capacity
  // (that is, if v.size() + NumberOfNewElements > v.capacity().
  // You can use v.reserve(NewCapacity) to ensure the vector has a capacity of at least NewCapacity.
  // Also remember that when the vector gets destroyed, the underlying array gets deleted as well.
  // None of the above conditions occurred during the period for which we use the conversion.

double OBForceFieldMMFF94GPU::EnergyElectrostatic() {

	double energy = 0.0;

	const int size = electrostaticCalculations.totalPairs();
	FloatType * ax = &electrostaticCalculations.ax[0];
	FloatType * ay = &electrostaticCalculations.ay[0];
	FloatType * az = &electrostaticCalculations.az[0];
	FloatType * bx = &electrostaticCalculations.bx[0];
	FloatType * by = &electrostaticCalculations.by[0];
	FloatType * bz = &electrostaticCalculations.bz[0];
	FloatType * qq = &electrostaticCalculations.qq[0];

	// calculate the bond vectors and determine its length (distance)
	FloatType abx, aby, abz;
	FloatType rab;

	for (int indx = 0; indx < size; indx++) {
		abx = ax[indx] - bx[indx];
		aby = ay[indx] - by[indx];
		abz = az[indx] - bz[indx];
		rab = sqrt(abx * abx + aby * aby + abz * abz);
		energy += qq[indx] / (rab + 0.05); // add 0.05 to avoid zero division
	} // for(electrostatic)

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



    bool OBForceFieldMMFF94GPU::SetupCalculations()
    {
        IF_OBFF_LOGLVL_LOW
          OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");

        bool gradients = false; // FIXME this parameter will come as an argument

    	bool setup;
    	setup = SetupBondCalculations(gradients);
    	setup &= SetupAngleAndStrBndCalculations(gradients);
    	setup &= SetupTorsionCalculations(gradients);
    	setup &= SetupOOPCalculations(gradients);
    	setup &= SetupVDWCalculations(gradients);
    	setup &= SetupElectrostaticCalculations(gradients);

    	return setup;
    }

    //
    // Bond Calculations
    //
    // no "step-down" procedure
    // MMFF part V - page 625 (empirical rule)
    //
    bool OBForceFieldMMFF94GPU::SetupBondCalculations(bool gradients)
    {
      OBFFParameter *parameter;
      OBAtom *a, *b, *c, *d;
      int type_a, type_b, type_c, type_d;
      bool found;
      int order;

      IF_OBFF_LOGLVL_LOW
        OBFFLog("SETTING UP BOND CALCULATIONS...\n");

      OBFFBondCalculationMMFF94 bondcalc;
      int bondtype;

      _bondcalculations.clear();

      FOR_BONDS_OF_MOL(bond, _mol) {
        a = bond->GetBeginAtom();
        b = bond->GetEndAtom();

        // skip this bond if the atoms are ignored
        if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
          continue;

        // if there are any groups specified, check if the two bond atoms are in a single intraGroup
        if (HasGroups()) {
          bool validBond = false;
          for (unsigned int i=0; i < _intraGroup.size(); ++i) {
            if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx())) {
              validBond = true;
              break;
            }
          }
          if (!validBond)
            continue;
        }

        bondtype = GetBondType(a, b);

        parameter = GetTypedParameter2Atom(bondtype, atoi(a->GetType()), atoi(b->GetType()), _ffbondparams); // from mmffbond.par
        if (parameter == NULL) {
          parameter = GetParameter2Atom(a->GetAtomicNum(), b->GetAtomicNum(), _ffbndkparams); // from mmffbndk.par - emperical rules
          if (parameter == NULL) {
            IF_OBFF_LOGLVL_LOW {
              // This should never happen
              snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR BOND %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
              OBFFLog(_logbuf);
            }
            return false;
          } else {
            IF_OBFF_LOGLVL_LOW {
              snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR BOND STRETCHING %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
              OBFFLog(_logbuf);
            }

            double rr, rr2, rr4, rr6;
            bondcalc.a = a;
            bondcalc.b = b;
            bondcalc.r0 = GetRuleBondLength(a, b);

            rr = parameter->_dpar[0] / bondcalc.r0; // parameter->_dpar[0]  = r0-ref
            rr2 = rr * rr;
            rr4 = rr2 * rr2;
            rr6 = rr4 * rr2;

            bondcalc.kb = parameter->_dpar[1] * rr6; // parameter->_dpar[1]  = kb-ref
            bondcalc.bt = bondtype;
            bondcalc.SetupPointers();

            _bondcalculations.push_back(bondcalc);
          }
        } else {
          bondcalc.a = a;
          bondcalc.b = b;
          bondcalc.kb = parameter->_dpar[0];
          bondcalc.r0 = parameter->_dpar[1];
          bondcalc.bt = bondtype;
          bondcalc.SetupPointers();

          _bondcalculations.push_back(bondcalc);
        }
      }

      return true;  // no proper return value !!
    }

      //
      // Angle Calculations
      //
      // MMFF part I - page 513 ("step-down" prodedure)
      // MMFF part I - page 519 (reference 68 is actually a footnote)
      // MMFF part V - page 627 (empirical rule)
      //
      // First try and find an exact match, if this fails, step down using the equivalences from mmffdef.par
      // five-stage protocol: 1-1-1, 2-2-2, 3-2-3, 4-2-4, 5-2-5
      // If this fails, use empirical rules
      // Since 1-1-1 = 2-2-2, we will only try 1-1-1 before going to 3-2-3
      //
      // Stretch-Bend Calculations
      //
    bool OBForceFieldMMFF94GPU::SetupAngleAndStrBndCalculations(bool gradients)
    {
      OBFFParameter *parameter;
      OBAtom *a, *b, *c, *d;
      int type_a, type_b, type_c, type_d;
      bool found;
      int order;

      IF_OBFF_LOGLVL_LOW
        OBFFLog("SETTING UP ANGLE & STRETCH-BEND CALCULATIONS...\n");

      OBFFAngleCalculationMMFF94 anglecalc;
      OBFFStrBndCalculationMMFF94 strbndcalc;
      int angletype, strbndtype, bondtype1, bondtype2;

      _anglecalculations.clear();
      _strbndcalculations.clear();

      FOR_ANGLES_OF_MOL(angle, _mol) {
        b = _mol.GetAtom((*angle)[0] + 1);
        a = _mol.GetAtom((*angle)[1] + 1);
        c = _mol.GetAtom((*angle)[2] + 1);

        type_a = atoi(a->GetType());
        type_b = atoi(b->GetType());
        type_c = atoi(c->GetType());

        // skip this angle if the atoms are ignored
        if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) || _constraints.IsIgnored(c->GetIdx()) )
          continue;

        // if there are any groups specified, check if the three angle atoms are in a single intraGroup
        if (HasGroups()) {
          bool validAngle = false;
          for (unsigned int i=0; i < _intraGroup.size(); ++i) {
            if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx()) &&
                _intraGroup[i].BitIsOn(c->GetIdx())) {
              validAngle = true;
              break;
            }
          }
          if (!validAngle)
            continue;
        }

        angletype = GetAngleType(a, b, c);
        strbndtype = GetStrBndType(a, b, c);
        bondtype1 = GetBondType(a, b);
        bondtype2 = GetBondType(b, c);

        if (HasLinSet(type_b)) {
          anglecalc.linear = true;
        } else {
          anglecalc.linear = false;
        }

        // try exact match
        parameter = GetTypedParameter3Atom(angletype, type_a, type_b, type_c, _ffangleparams);
        if (parameter == NULL) // try 3-2-3
          parameter = GetTypedParameter3Atom(angletype, EqLvl3(type_a), type_b, EqLvl3(type_c), _ffangleparams);
        if (parameter == NULL) // try 4-2-4
          parameter = GetTypedParameter3Atom(angletype, EqLvl4(type_a), type_b, EqLvl4(type_c), _ffangleparams);
        if (parameter == NULL) // try 5-2-5
          parameter = GetTypedParameter3Atom(angletype, EqLvl5(type_a), type_b, EqLvl5(type_c), _ffangleparams);

        if (parameter) {
          anglecalc.ka = parameter->_dpar[0];
          anglecalc.theta0 = parameter->_dpar[1];
          strbndcalc.theta0 = parameter->_dpar[1]; // **
        } else {
          IF_OBFF_LOGLVL_LOW {
            snprintf(_logbuf, BUFF_SIZE, "   USING DEFAULT ANGLE FOR %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
            snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR ANGLE BENDING %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
            OBFFLog(_logbuf);
          }

          anglecalc.ka = 0.0;
          anglecalc.theta0 = 120.0;

          if (GetCrd(type_b) == 4)
            anglecalc.theta0 = 109.45;

          if ((GetCrd(type_b) == 2) && b->IsOxygen())
            anglecalc.theta0 = 105.0;

          if (b->GetAtomicNum() > 10)
            anglecalc.theta0 = 95.0;

          if (HasLinSet(type_b))
            anglecalc.theta0 = 180.0;

          if ((GetCrd(type_b) == 3) && (GetVal(type_b) == 3) && !GetMltb(type_b)) {
            if (b->IsNitrogen()) {
              anglecalc.theta0 = 107.0;
            } else {
              anglecalc.theta0 = 92.0;
            }
          }

          if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
            anglecalc.theta0 = 60.0;

          if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
            anglecalc.theta0 = 90.0;

          strbndcalc.theta0 = anglecalc.theta0; // **
        }

        // empirical rule for 0-b-0 and standard angles
        if (anglecalc.ka == 0.0) {
          IF_OBFF_LOGLVL_LOW {
            snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR ANGLE BENDING FORCE CONSTANT %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
            OBFFLog(_logbuf);
          }

          double beta, Za, Zc, Cb, r0ab, r0bc, theta, theta2, D, rr, rr2;
          Za = GetZParam(a);
          Cb = GetCParam(b); // Fixed typo -- PR#2741658
          Zc = GetZParam(c);

          r0ab = GetBondLength(a, b);
          r0bc = GetBondLength(b, c);
          rr = r0ab + r0bc;
          rr2 = rr * rr;
          D = (r0ab - r0bc) / rr2;

          theta = anglecalc.theta0;
          theta2 = theta * theta;

          beta = 1.75;
          if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
            beta = 0.85 * beta;
          if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
            beta = 0.05 * beta;

          // Theta2 is in Degrees^2, but parameters are expecting radians
          // PR#2741669
          anglecalc.ka = (beta * Za * Cb * Zc * exp(-2 * D)) / (rr * theta2 * DEG_TO_RAD * DEG_TO_RAD);
        }

        anglecalc.a = a;
        anglecalc.b = b;
        anglecalc.c = c;
        anglecalc.at = angletype;

        anglecalc.SetupPointers();
        _anglecalculations.push_back(anglecalc);

        if (anglecalc.linear)
          continue;

        parameter = GetTypedParameter3Atom(strbndtype, type_a, type_b, type_c, _ffstrbndparams);
        if (parameter == NULL) {
          int rowa, rowb, rowc;

          rowa = GetElementRow(a);
          rowb = GetElementRow(b);
          rowc = GetElementRow(c);

          parameter = GetParameter3Atom(rowa, rowb, rowc, _ffdfsbparams);

          if (parameter == NULL) {
            // This should never happen
            IF_OBFF_LOGLVL_LOW {
              snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR STRETCH-BEND %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
              OBFFLog(_logbuf);
            }
            return false;
          }

          if (rowa == parameter->a) {
            strbndcalc.kbaABC = parameter->_dpar[0];
            strbndcalc.kbaCBA = parameter->_dpar[1];
          } else {
            strbndcalc.kbaABC = parameter->_dpar[1];
            strbndcalc.kbaCBA = parameter->_dpar[0];
          }
        } else {
          if (type_a == parameter->a) {
            strbndcalc.kbaABC = parameter->_dpar[0];
            strbndcalc.kbaCBA = parameter->_dpar[1];
          } else {
            strbndcalc.kbaABC = parameter->_dpar[1];
            strbndcalc.kbaCBA = parameter->_dpar[0];
          }
        }

        strbndcalc.rab0 = GetBondLength(a, b);
        strbndcalc.rbc0 = GetBondLength(b ,c);
        strbndcalc.a = a;
        strbndcalc.b = b;
        strbndcalc.c = c;
        strbndcalc.sbt = strbndtype;
        strbndcalc.SetupPointers();

        _strbndcalculations.push_back(strbndcalc);

      }

      return true; //FIXME here should return something meaningful or declare method as void
    }

      //
      // Torsion Calculations
      //
      // MMFF part I - page 513 ("step-down" prodedure)
      // MMFF part I - page 519 (reference 68 is actually a footnote)
      // MMFF part IV - page 631 (empirical rule)
      //
      // First try and find an exact match, if this fails, step down using the equivalences from mmffdef.par
      // five-stage protocol: 1-1-1-1, 2-2-2-2, 3-2-2-5, 5-2-2-3, 5-2-2-5
      // If this fails, use empirical rules
      // Since 1-1-1-1 = 2-2-2-2, we will only try 1-1-1-1 before going to 3-2-2-5
      //
    bool OBForceFieldMMFF94GPU::SetupTorsionCalculations(bool gradients)
    {
      OBFFParameter *parameter;
      OBAtom *a, *b, *c, *d;
      int type_a, type_b, type_c, type_d;
      bool found;
      int order;

      IF_OBFF_LOGLVL_LOW
        OBFFLog("SETTING UP TORSION CALCULATIONS...\n");

      OBFFTorsionCalculationMMFF94 torsioncalc;
      int torsiontype;

      _torsioncalculations.clear();

      FOR_TORSIONS_OF_MOL(t, _mol) {
        a = _mol.GetAtom((*t)[0] + 1);
        b = _mol.GetAtom((*t)[1] + 1);
        c = _mol.GetAtom((*t)[2] + 1);
        d = _mol.GetAtom((*t)[3] + 1);

        type_a = atoi(a->GetType());
        type_b = atoi(b->GetType());
        type_c = atoi(c->GetType());
        type_d = atoi(d->GetType());

        // skip this torsion if the atoms are ignored
        if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
             _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) )
          continue;

        // if there are any groups specified, check if the four torsion atoms are in a single intraGroup
        if (HasGroups()) {
          bool validTorsion = false;
          for (unsigned int i=0; i < _intraGroup.size(); ++i) {
            if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx()) &&
                _intraGroup[i].BitIsOn(c->GetIdx()) && _intraGroup[i].BitIsOn(d->GetIdx())) {
              validTorsion = true;
              break;
            }
          }
          if (!validTorsion)
            continue;
        }

        torsiontype = GetTorsionType(a, b, c, d);
        // CXT = MC*(J*MA**3 + K*MA**2 + I*MA + L) + TTijkl  MC = 6, MA = 136
        order = (type_c*2515456 + type_b*18496 + type_d*136 + type_a)
          - (type_b*2515456 + type_c*18496 + type_a*136 + type_d);

        if (order >= 0) {
          // try exact match
          parameter = GetTypedParameter4Atom(torsiontype, type_a, type_b, type_c, type_d, _fftorsionparams);
          if (parameter == NULL) // try 3-2-2-5
            parameter = GetTypedParameter4Atom(torsiontype, EqLvl3(type_a), type_b, type_c, EqLvl5(type_d), _fftorsionparams);
          if (parameter == NULL) // try 5-2-2-3
            parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_a), type_b, type_c, EqLvl3(type_d), _fftorsionparams);
          if (parameter == NULL) // try 5-2-2-5
            parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_a), type_b, type_c, EqLvl5(type_d), _fftorsionparams);
        } else {
          // try exact match
          parameter = GetTypedParameter4Atom(torsiontype, type_d, type_c, type_b, type_a, _fftorsionparams);
          if (parameter == NULL) // try 3-2-2-5
            parameter = GetTypedParameter4Atom(torsiontype, EqLvl3(type_d), type_c, type_b, EqLvl5(type_a), _fftorsionparams);
          if (parameter == NULL) // try 5-2-2-3
            parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_d), type_c, type_b, EqLvl3(type_a), _fftorsionparams);
          if (parameter == NULL) // try 5-2-2-5
            parameter = GetTypedParameter4Atom(torsiontype, EqLvl5(type_d), type_c, type_b, EqLvl5(type_a), _fftorsionparams);
        }

        if (parameter) {
          torsioncalc.v1 = parameter->_dpar[0];
          torsioncalc.v2 = parameter->_dpar[1];
          torsioncalc.v3 = parameter->_dpar[2];
        } else {
          bool found_rule = false;

          //IF_OBFF_LOGLVL_LOW {
          //  snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR TORSION FORCE CONSTANT %d-%d-%d-%d (IDX)...\n",
          //    a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());
          //  OBFFLog(_logbuf);
          //}

          // rule (a) page 631
          if (HasLinSet(type_b) || HasLinSet(type_c))
            continue;

          // rule (b) page 631
          if (b->GetBond(c)->IsAromatic()) {
            double Ub, Uc, pi_bc, beta;
            Ub = GetUParam(b);
            Uc = GetUParam(c);

            if (!HasPilpSet(type_b) && !HasPilpSet(type_c))
              pi_bc = 0.5;
            else
              pi_bc = 0.3;

            if (((GetVal(type_b) == 3) && (GetVal(type_c) == 4)) ||
                ((GetVal(type_b) == 4) && (GetVal(type_c) == 3)))
              beta = 3.0;
            else
              beta = 6.0;

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
            torsioncalc.v3 = 0.0;
            found_rule = true;
          } else {
            // rule (c) page 631
         	  double Ub, Uc, pi_bc, beta;
            Ub = GetUParam(b);
            Uc = GetUParam(c);

            if (((GetMltb(type_b) == 2) && (GetMltb(type_c) == 2)) && a->GetBond(b)->IsDouble())
              pi_bc = 1.0;
            else
              pi_bc = 0.4;

            beta = 6.0;
            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
            torsioncalc.v3 = 0.0;
            found_rule = true;
          }

          // rule (d) page 632
          if (!found_rule)
            if (((GetCrd(type_b) == 4) && (GetCrd(type_c) == 4))) {
              double Vb, Vc;
              Vb = GetVParam(b);
              Vc = GetVParam(c);

              torsioncalc.v1 = 0.0;
              torsioncalc.v2 = 0.0;
              torsioncalc.v3 = sqrt(Vb * Vc) / 9.0;
              found_rule = true;
            }

          // rule (e) page 632
          if (!found_rule)
            if (((GetCrd(type_b) == 4) && (GetCrd(type_c) != 4))) {
              if (GetCrd(type_c) == 3) // case (1)
                if ((GetVal(type_c) == 4) || (GetVal(type_c) == 34) || (GetMltb(type_c) != 0))
                  continue;

              if (GetCrd(type_c) == 2) // case (2)
                if ((GetVal(type_c) == 3) || (GetMltb(type_c) != 0))
                  continue;

              // case (3) saturated bonds -- see rule (h)
            }

          // rule (f) page 632
          if (!found_rule)
            if (((GetCrd(type_b) != 4) && (GetCrd(type_c) == 4))) {
              if (GetCrd(type_b) == 3) // case (1)
                if ((GetVal(type_b) == 4) || (GetVal(type_b) == 34) || (GetMltb(type_b) != 0))
                  continue;

              if (GetCrd(type_b) == 2) // case (2)
                if ((GetVal(type_b) == 3) || (GetMltb(type_b) != 0))
                  continue;

              // case (3) saturated bonds
            }

          // rule (g) page 632
          if (!found_rule)
            if (b->GetBond(c)->IsSingle() && (
                                              (GetMltb(type_b) && GetMltb(type_c)) ||
                                              (GetMltb(type_b) && HasPilpSet(type_c)) ||
                                              (GetMltb(type_c) && HasPilpSet(type_b))  )) {
              if (HasPilpSet(type_b) && HasPilpSet(type_c)) // case (1)
                continue;

              double Ub, Uc, pi_bc, beta;
              Ub = GetUParam(b);
              Uc = GetUParam(c);
              beta = 6.0;

              if (HasPilpSet(type_b) && GetMltb(type_c)) { // case (2)
                if (GetMltb(type_c) == 1)
                  pi_bc = 0.5;
                else if ((GetElementRow(b) == 1) && (GetElementRow(c) == 1))
                  pi_bc = 0.3;
                else
                  pi_bc = 0.15;
                found_rule = true;
              }

              if (HasPilpSet(type_c) && GetMltb(type_b)) { // case (3)
                if (GetMltb(type_b) == 1)
                  pi_bc = 0.5;
                else if ((GetElementRow(b) == 1) && (GetElementRow(c) == 1))
                  pi_bc = 0.3;
                else
                  pi_bc = 0.15;
                found_rule = true;
              }

              if (!found_rule)
                if (((GetMltb(type_b) == 1) || (GetMltb(type_c) == 1)) && (!b->IsCarbon() || !c->IsCarbon())) {
                  pi_bc = 0.4;
                  found_rule = true;
                }

              if (!found_rule)
                pi_bc = 0.15;

              torsioncalc.v1 = 0.0;
              torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
              torsioncalc.v3 = 0.0;
              found_rule = true;
            }

          // rule (h) page 632
          if (!found_rule) {
            if ((b->IsOxygen() || b->IsSulfur()) && (c->IsOxygen() || c->IsSulfur())) {
              double Wb, Wc;

              if (b->IsOxygen()) {
                Wb = 2.0;
              }
              else {
                Wb = 8.0;
              }

              if (c->IsOxygen()) {
                Wc = 2.0;
              }
              else {
                Wc = 8.0;
              }

              torsioncalc.v1 = 0.0;
              torsioncalc.v2 = -sqrt(Wb * Wc);
              torsioncalc.v3 = 0.0;
            } else {
              double Vb, Vc, Nbc;
              Vb = GetVParam(b);
              Vc = GetVParam(c);

              IF_OBFF_LOGLVL_LOW {
                snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR TORSION FORCE CONSTANT %d-%d-%d-%d (IDX)...\n",
                        a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());
                OBFFLog(_logbuf);
              }

              Nbc = GetCrd(type_b) * GetCrd(type_c);

              torsioncalc.v1 = 0.0;
              torsioncalc.v2 = 0.0;
              torsioncalc.v3 = sqrt(Vb * Vc) / Nbc;
            }
          }
        }

        torsioncalc.a = a;
        torsioncalc.b = b;
        torsioncalc.c = c;
        torsioncalc.d = d;
        torsioncalc.SetupPointers();
        torsioncalc.tt = torsiontype;

        _torsioncalculations.push_back(torsioncalc);
      }

      return true;
    }

      //
      // Out-Of-Plane Calculations
      //
    bool OBForceFieldMMFF94GPU::SetupOOPCalculations(bool gradients)
    {
      OBFFParameter *parameter;
      OBAtom *a, *b, *c, *d;
      int type_a, type_b, type_c, type_d;
      bool found;
      int order;

    IF_OBFF_LOGLVL_LOW
        OBFFLog("SETTING UP OOP CALCULATIONS...\n");

      OBFFOOPCalculationMMFF94 oopcalc;

      _oopcalculations.clear();

      FOR_ATOMS_OF_MOL(atom, _mol) {
        b = (OBAtom*) &*atom;

        found = false;

        type_b = atoi(b->GetType());

        for (unsigned int idx=0; idx < _ffoopparams.size(); idx++) {
          if (type_b == _ffoopparams[idx].b) {
            a = NULL;
            c = NULL;
            d = NULL;

            FOR_NBORS_OF_ATOM(nbr, b) {
              if (a ==NULL)
                a = (OBAtom*) &*nbr;
              else if (c == NULL)
                c = (OBAtom*) &*nbr;
              else
                d = (OBAtom*) &*nbr;
            }

            if ((a == NULL) || (c == NULL) || (d == NULL))
              break;

            type_a = atoi(a->GetType());
            type_c = atoi(c->GetType());
            type_d = atoi(d->GetType());

            // skip this oop if the atoms are ignored
            if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
                 _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) )
              continue;

            // if there are any groups specified, check if the four oop atoms are in a single intraGroup
            if (HasGroups()) {
              bool validOOP = false;
              for (unsigned int i=0; i < _intraGroup.size(); ++i) {
                if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx()) &&
                    _intraGroup[i].BitIsOn(c->GetIdx()) && _intraGroup[i].BitIsOn(d->GetIdx())) {
                  validOOP = true;
                  break;
                }
              }
              if (!validOOP)
                continue;
            }

            if (((type_a == _ffoopparams[idx].a) && (type_c == _ffoopparams[idx].c) && (type_d == _ffoopparams[idx].d)) ||
                ((type_c == _ffoopparams[idx].a) && (type_a == _ffoopparams[idx].c) && (type_d == _ffoopparams[idx].d)) ||
                ((type_c == _ffoopparams[idx].a) && (type_d == _ffoopparams[idx].c) && (type_a == _ffoopparams[idx].d)) ||
                ((type_d == _ffoopparams[idx].a) && (type_c == _ffoopparams[idx].c) && (type_a == _ffoopparams[idx].d)) ||
                ((type_a == _ffoopparams[idx].a) && (type_d == _ffoopparams[idx].c) && (type_c == _ffoopparams[idx].d)) ||
                ((type_d == _ffoopparams[idx].a) && (type_a == _ffoopparams[idx].c) && (type_c == _ffoopparams[idx].d)))
              {
                found = true;

                oopcalc.koop = _ffoopparams[idx]._dpar[0];

                // A-B-CD || C-B-AD  PLANE = ABC
                oopcalc.a = a;
                oopcalc.b = b;
                oopcalc.c = c;
                oopcalc.d = d;

                oopcalc.SetupPointers();
                _oopcalculations.push_back(oopcalc);

                // C-B-DA || D-B-CA  PLANE BCD
                oopcalc.a = d;
                oopcalc.d = a;

                oopcalc.SetupPointers();
                _oopcalculations.push_back(oopcalc);

                // A-B-DC || D-B-AC  PLANE ABD
                oopcalc.a = a;
                oopcalc.c = d;
                oopcalc.d = c;

                oopcalc.SetupPointers();
                _oopcalculations.push_back(oopcalc);
              }

            if ((_ffoopparams[idx].a == 0) && (_ffoopparams[idx].c == 0) && (_ffoopparams[idx].d == 0) && !found) // *-XX-*-*
              {
                oopcalc.koop = _ffoopparams[idx]._dpar[0];

                // A-B-CD || C-B-AD  PLANE = ABC
                oopcalc.a = a;
                oopcalc.b = b;
                oopcalc.c = c;
                oopcalc.d = d;

                oopcalc.SetupPointers();
                _oopcalculations.push_back(oopcalc);

                // C-B-DA || D-B-CA  PLANE BCD
                oopcalc.a = d;
                oopcalc.d = a;

                oopcalc.SetupPointers();
                _oopcalculations.push_back(oopcalc);

                // A-B-DC || D-B-AC  PLANE ABD
                oopcalc.a = a;
                oopcalc.c = d;
                oopcalc.d = c;

                oopcalc.SetupPointers();
                _oopcalculations.push_back(oopcalc);
              }
          }
        }
      }

      return true; //return something meaningful or change to void
    }

      //
      // VDW Calculations
      //
    bool OBForceFieldMMFF94GPU::SetupVDWCalculations(bool gradients)
    {
      OBFFParameter *parameter;
      OBAtom *a, *b, *c, *d;
      int type_a, type_b, type_c, type_d;
      bool found;
      int order;


      IF_OBFF_LOGLVL_LOW
        OBFFLog("SETTING UP VAN DER WAALS CALCULATIONS...\n");

      vdwCalculations.reset();

      int pairIndex = -1;
      double *pos_a, *pos_b;
      FOR_PAIRS_OF_MOL(p, _mol) {
        ++pairIndex;
        a = _mol.GetAtom((*p)[0]);
        b = _mol.GetAtom((*p)[1]);

        // skip this vdw if the atoms are ignored
        if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
          continue;

        // if there are any groups specified, check if the two atoms are in a single _interGroup or if
        // two two atoms are in one of the _interGroups pairs.
        if (HasGroups()) {
          bool validVDW = false;
          for (unsigned int i=0; i < _interGroup.size(); ++i) {
            if (_interGroup[i].BitIsOn(a->GetIdx()) && _interGroup[i].BitIsOn(b->GetIdx())) {
              validVDW = true;
              break;
            }
          }
          if (!validVDW) {
            for (unsigned int i=0; i < _interGroups.size(); ++i) {
              if (_interGroups[i].first.BitIsOn(a->GetIdx()) && _interGroups[i].second.BitIsOn(b->GetIdx())) {
                validVDW = true;
                break;
              }
              if (_interGroups[i].first.BitIsOn(b->GetIdx()) && _interGroups[i].second.BitIsOn(a->GetIdx())) {
                validVDW = true;
                break;
              }
            }
          }

          if (!validVDW)
            continue;
        }

        OBFFParameter *parameter_a, *parameter_b;
        parameter_a = GetParameter1Atom(atoi(a->GetType()), _ffvdwparams);
        parameter_b = GetParameter1Atom(atoi(b->GetType()), _ffvdwparams);
        if ((parameter_a == NULL) || (parameter_b == NULL)) {
          IF_OBFF_LOGLVL_LOW {
            snprintf(_logbuf, BUFF_SIZE, "   COULD NOT FIND VAN DER WAALS PARAMETERS FOR %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
            OBFFLog(_logbuf);
          }

          return false;
        }

        int aDA, bDA; // hydrogen donor/acceptor (A=1, D=2, neither=0)
        double rab, epsilon, alpha_a, alpha_b, Na, Nb, Aa, Ab, Ga, Gb;
        double R_AB;

        alpha_a = parameter_a->_dpar[0];
        Na = parameter_a->_dpar[1];
        Aa = parameter_a->_dpar[2];
        Ga = parameter_a->_dpar[3];
        aDA = parameter_a->_ipar[0];

        alpha_b = parameter_b->_dpar[0];
        Nb = parameter_b->_dpar[1];
        Ab = parameter_b->_dpar[2];
        Gb = parameter_b->_dpar[3];
        bDA = parameter_b->_ipar[0];

        //these calculations only need to be done once for each pair,
        //we do them now and save them for later use
        double R_AA, R_BB, g_AB, g_AB2;
        double sqrt_a, sqrt_b;

        R_AA = Aa * pow(alpha_a, 0.25);
        R_BB = Ab * pow(alpha_b, 0.25);
        sqrt_a = sqrt(alpha_a / Na);
        sqrt_b = sqrt(alpha_b / Nb);

        if (aDA == 1) { // hydrogen bond donor
          R_AB = 0.5 * (R_AA + R_BB);

          if (bDA == 2) { // hydrogen bond acceptor
            epsilon = 0.5 * (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / pow(R_AB,6));
            // R_AB is scaled to 0.8 for D-A interactions.
            // NOTE!! The value used in the calculation of epsilon is not scaled.
            // R_AB7 however uses the new scaled value of R_AB
            R_AB = 0.8 * R_AB;
          } else {
            epsilon = (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / pow(R_AB,6));
          }

        } else if (bDA == 1) { // hydrogen bond donor
          R_AB = 0.5 * (R_AA + R_BB);

          if (aDA == 2) { // hydrogen bond acceptor
            epsilon = 0.5 * (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / pow(R_AB,6));
            // R_AB is scaled to 0.8 for D-A interactions.
            // NOTE!! The value used in the calculation of epsilon is not scaled.
            // R_AB7 however uses the new scaled value of R_AB
            R_AB = 0.8 * R_AB;
          } else {
            epsilon = (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / pow(R_AB,6));
          }

        } else {
          g_AB = (R_AA - R_BB) / ( R_AA + R_BB);
          g_AB2 = g_AB * g_AB;
          R_AB =  0.5 * (R_AA + R_BB) * (1.0 + 0.2 * (1.0 - exp(-12.0 * g_AB2)));
          epsilon = (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / pow(R_AB,6));
        }

		pos_a = a->GetCoordinate();
		pos_b = b->GetCoordinate();

		vdwCalculations.ax.push_back(pos_a[0]);
		vdwCalculations.ay.push_back(pos_a[1]);
		vdwCalculations.az.push_back(pos_a[2]);
		vdwCalculations.bx.push_back(pos_b[0]);
		vdwCalculations.by.push_back(pos_b[1]);
		vdwCalculations.bz.push_back(pos_b[2]);
		vdwCalculations.RAB.push_back(R_AB);
		vdwCalculations.epsilon.push_back(epsilon);
		vdwCalculations.indexA.push_back(a->GetIdx());
		vdwCalculations.indexB.push_back(b->GetIdx());
		vdwCalculations.pairIndex.push_back(pairIndex);

      }

      return true; // return something meaningful or change to void
    }

      //
      // Electrostatic Calculations
      //
    bool OBForceFieldMMFF94GPU::SetupElectrostaticCalculations(bool gradients)
    {
      OBFFParameter *parameter;
      OBAtom *a, *b, *c, *d;
      int type_a, type_b, type_c, type_d;
      bool found;
      int order;


      IF_OBFF_LOGLVL_LOW
        OBFFLog("SETTING UP ELECTROSTATIC CALCULATIONS...\n");

      electrostaticCalculations.reset();

      int pairIndex = -1;
      double qq = 0;
      double *pos_a, *pos_b;

      FOR_PAIRS_OF_MOL(p, _mol) {
        ++pairIndex;
        a = _mol.GetAtom((*p)[0]);
        b = _mol.GetAtom((*p)[1]);

        // skip this ele if the atoms are ignored
        if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
          continue;

        // if there are any groups specified, check if the two atoms are in a single _interGroup or if
        // two two atoms are in one of the _interGroups pairs.
        if (HasGroups()) {
          bool validEle = false;
          for (unsigned int i=0; i < _interGroup.size(); ++i) {
            if (_interGroup[i].BitIsOn(a->GetIdx()) && _interGroup[i].BitIsOn(b->GetIdx())) {
              validEle = true;
              break;
            }
          }
          if (!validEle) {
            for (unsigned int i=0; i < _interGroups.size(); ++i) {
              if (_interGroups[i].first.BitIsOn(a->GetIdx()) && _interGroups[i].second.BitIsOn(b->GetIdx())) {
                validEle = true;
                break;
              }
              if (_interGroups[i].first.BitIsOn(b->GetIdx()) && _interGroups[i].second.BitIsOn(a->GetIdx())) {
                validEle = true;
                break;
              }
            }
          }

          if (!validEle)
            continue;
        }

        qq = 332.0716 * a->GetPartialCharge() * b->GetPartialCharge();

		if (qq) {
			pos_a = a->GetCoordinate();
			pos_b = b->GetCoordinate();

			// 1-4 scaling
			if (a->IsOneFour(b))
				qq *= 0.75;

			electrostaticCalculations.ax.push_back(pos_a[0]);
			electrostaticCalculations.ay.push_back(pos_a[1]);
			electrostaticCalculations.az.push_back(pos_a[2]);
			electrostaticCalculations.bx.push_back(pos_b[0]);
			electrostaticCalculations.by.push_back(pos_b[1]);
			electrostaticCalculations.bz.push_back(pos_b[2]);
			electrostaticCalculations.qq.push_back(qq);
			electrostaticCalculations.indexA.push_back(a->GetIdx());
			electrostaticCalculations.indexB.push_back(b->GetIdx());
			electrostaticCalculations.pairIndex.push_back(pairIndex);

		} //if(qq)
      }//for(pairs)

      return true;
    }

} // end namespace OpenBabel

//! \file forcefieldmmff94gpu.cpp
//! \brief MMFF94 force field

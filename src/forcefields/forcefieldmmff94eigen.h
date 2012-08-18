/**********************************************************************
forcefieldmmff94eigen.h - MMFF94 using Eigen

Based on forcefieldmmff94.h
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

Modifications to run using Eigen acceleration
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

#include <vector>
#include <string>
#include <map>

#include <openbabel/parsmart.h>
#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

#include <Eigen/Dense>

#include "forcefieldmmff94.h"

#include <openbabel/benchmark-utils.hpp>
//#include <openbabel/mapkeys.h>


using namespace std; // in order to be able to succesfully use cout

namespace OpenBabel
{

class OBFFTorsionCalculationMMFF94Eigen {
public:
	std::vector<Eigen::Vector3d> posAVector;
	std::vector<Eigen::Vector3d> posBVector;
	std::vector<Eigen::Vector3d> posCVector;
	std::vector<Eigen::Vector3d> posDVector;
	std::vector<Eigen::Vector3d> velocityVector;
    std::vector<int> torsionTypeVector; //torsiontype (TTIJKL)


	void reset(){
		posAVector.clear();
		posBVector.clear();
		posCVector.clear();
		posDVector.clear();
		velocityVector.clear();
		torsionTypeVector.clear();
	}//reset()

	int totalCalcs(){
		return torsionTypeVector.size(); // all vectors in torsion calculations are same size
	}
};


class OBFFOOPCalculationMMFF94Eigen {
public:
	std::vector<Eigen::Vector3d> posAVector;
	std::vector<Eigen::Vector3d> posBVector;
	std::vector<Eigen::Vector3d> posCVector;
	std::vector<Eigen::Vector3d> posDVector;
	std::vector<double> koopVector;

	void reset(){
		posAVector.clear();
		posBVector.clear();
		posCVector.clear();
		posDVector.clear();
		koopVector.clear();
	}//reset()

	int totalCalcs(){
		return koopVector.size(); // all vectors in OOP are same size
	}
};


class OBFFElectrostaticCalculationMMFF94Eigen {
public:
	std::vector<Eigen::Vector3d> posAVector;
	std::vector<Eigen::Vector3d> posBVector;
	std::vector<double> qqVector;
	std::vector<int> indexA; // index atom A
	std::vector<int> indexB; // index atom B
	std::vector<int> pairIndex; // index into iteration using FOR_PAIRS_OF_MOL(..., _mol)

	void reset(){
		posAVector.clear();
		posBVector.clear();
		qqVector.clear();
		indexA.clear();
		indexB.clear();
		pairIndex.clear();
	}//reset()

	int totalPairs(){
		return pairIndex.size();
	}
};

class OBFFVDWCalculationMMFF94Eigen {
public:
    std::vector<Eigen::Vector3d> posAVector;
    std::vector<Eigen::Vector3d> posBVector;
    std::vector<double> RABVector;
    std::vector<double> RAB7Vector;
    std::vector<double> epsilonVector;
	std::vector<int> indexA; // index atom A
	std::vector<int> indexB; // index atom B
	std::vector<int> pairIndex; // index into iteration using FOR_PAIRS_OF_MOL(..., _mol)

	void reset(){
		posAVector.clear();
		posBVector.clear();
		RABVector.clear();
		RAB7Vector.clear();
		epsilonVector.clear();
		indexA.clear();
		indexB.clear();
		pairIndex.clear();
	}//reset()

	int totalPairs(){
		return pairIndex.size();
	}

};



  // Class OBForceFieldMMFF94Eigen
  // class introduction in forcefieldmmff94eigen.cpp
  class OBForceFieldMMFF94Eigen: public OBForceFieldMMFF94
  {
  protected:

      // OBFFXXXCalculationYYY structures to contain the calculations
//      OBFFBondCalculationMMFF94Eigen bondCalculations;
//      OBFFAngleCalculationMMFF94Eigen angleCalculations;
//      OBFFStrBndCalculationMMFF94Eigen strbndCalculations;
      OBFFTorsionCalculationMMFF94Eigen torsionCalculations;
      OBFFOOPCalculationMMFF94Eigen oopCalculations;
	  OBFFVDWCalculationMMFF94Eigen vdwCalculations;
	  OBFFElectrostaticCalculationMMFF94Eigen electrostaticCalculations;

    public:
      //! Constructor
      explicit OBForceFieldMMFF94Eigen(const char* ID, bool IsDefault=true) : OBForceFieldMMFF94(ID, IsDefault)
      {
    	  // The initialization of parameters default values is carried out in the parent class OBForceFieldMMFF94

      }

      //! Destructor
      virtual ~OBForceFieldMMFF94Eigen();

      //! Assignment
      OBForceFieldMMFF94Eigen &operator = (OBForceFieldMMFF94Eigen &);

      //! Returns total energy
      double Energy(bool gradients = true);
      //! Returns the bond stretching energy
      double EnergyBond();
      //! Returns the angle bending energy
      double EnergyAngle();
      //! Returns the stretch-bend energy
      double EnergyStrBnd();
      //! Returns the torsional energy
      double EnergyTorsion();
      //! Returns the out-of-plane bending energy
      double EnergyOOP();
      //! Returns the Van der Waals energy (Buckingham potential)
      double EnergyVDW();
      //! Returns the dipole-dipole interaction energy
      double EnergyElectrostatic();

      //! Setup the model calculation parameters
      bool SetupCalculations();
//      bool SetupCalculations(bool gradients = false);

      //! Setup the bond stretching energy calculations
      bool SetupBondCalculations(bool gradients);
      //! Setup the angle bending energy and the stretch-bend energy calculations
      bool SetupAngleAndStrBndCalculations(bool gradients);
      //! Setup the torsional energy calculations
      bool SetupTorsionCalculations(bool gradients);
      //! Setup the out-of-plane bending energy calculations
      bool SetupOOPCalculations(bool gradients);
      //! Setup the Van der Waals energy (Buckingham potential) calculations
      bool SetupVDWCalculations(bool gradients);
      //! Setup the dipole-dipole interaction energy calculations
      bool SetupElectrostaticCalculations(bool gradients);

  }; // class OBForceFieldMMFF94Eigen

}// namespace OpenBabel

//! \file forcefieldmmff94eigen.h
//! \brief MMFF94 force field


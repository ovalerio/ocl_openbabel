/**********************************************************************
forcefieldmmff94gpu.h - MMFF94 using GPU

Based on forcefieldmmff94.h
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>

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

#include <vector>
#include <string>
#include <map>

#include <openbabel/parsmart.h>
#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

#include "forcefieldmmff94.h"

#define FloatType double

using namespace std; // in order to be able to succesfully use cout

namespace OpenBabel
{

/*  compiler has better idea of doing vectorization when using Structure of Arrays */
class OBFFElectrostaticCalculationMMFF94GPU {
public:
	std::vector<FloatType> ax;
	std::vector<FloatType> ay;
	std::vector<FloatType> az;
	std::vector<FloatType> bx;
	std::vector<FloatType> by;
	std::vector<FloatType> bz;
	std::vector<FloatType> qq;
	std::vector<int> indexA; // index atom A
	std::vector<int> indexB; // index atom B
	std::vector<int> pairIndex; // index into iteration using FOR_PAIRS_OF_MOL(..., _mol)

//	//! Constructor
//	explicit OBFFElectrostaticCalculationMMFF94GPU() {
//
//	}
//
//	//! Destructor
//	virtual ~OBFFElectrostaticCalculationMMFF94GPU();

//		template<bool> void Compute();

	void reset(){
		ax.clear();
		ay.clear();
		az.clear();
		bx.clear();
		by.clear();
		bz.clear();
		qq.clear();
		indexA.clear();
		indexB.clear();
		pairIndex.clear();
	}//reset()

	int totalPairs(){
		return pairIndex.size();
	}
};

class OBFFVDWCalculationMMFF94GPU {
public:
	std::vector<FloatType> ax;
	std::vector<FloatType> ay;
	std::vector<FloatType> az;
	std::vector<FloatType> bx;
	std::vector<FloatType> by;
	std::vector<FloatType> bz;
	std::vector<FloatType> RAB;
	std::vector<FloatType> epsilon;
	std::vector<int> indexA; // index atom A
	std::vector<int> indexB; // index atom B
	std::vector<int> pairIndex; // index into iteration using FOR_PAIRS_OF_MOL(..., _mol)

	//! Constructor
//	explicit OBFFVDWCalculationMMFF94GPU() {
//
//	}
//
//	//! Destructor
//	virtual ~OBFFVDWCalculationMMFF94GPU();

//		template<bool> void Compute();

	void reset(){
		ax.clear();
		ay.clear();
		az.clear();
		bx.clear();
		by.clear();
		bz.clear();
		RAB.clear();
		epsilon.clear();
		indexA.clear();
		indexB.clear();
		pairIndex.clear();
	}//reset()

	int totalPairs(){
		return pairIndex.size();
	}

};


  // Class OBForceFieldMMFF94GPU
  // class introduction in forcefieldmmff94gpu.cpp
  class OBForceFieldMMFF94GPU: public OBForceFieldMMFF94
  {
  protected:

      // OBFFXXXCalculationYYY structures to contain the calculations
//      std::vector<OBFFBondCalculationMMFF94>          _bondcalculations;
//      std::vector<OBFFAngleCalculationMMFF94>         _anglecalculations;
//      std::vector<OBFFStrBndCalculationMMFF94>        _strbndcalculations;
//      std::vector<OBFFTorsionCalculationMMFF94>       _torsioncalculations;
//      std::vector<OBFFOOPCalculationMMFF94>           _oopcalculations;
	  OBFFVDWCalculationMMFF94GPU vdwCalculations;
	  OBFFElectrostaticCalculationMMFF94GPU electrostaticCalculations;

    public:
      //! Constructor
      explicit OBForceFieldMMFF94GPU(const char* ID, bool IsDefault=true) : OBForceFieldMMFF94(ID, IsDefault)
      {
    	  // The initialization of parameters default values is carried out in the parent class OBForceFieldMMFF94

      }

      //! Destructor
      virtual ~OBForceFieldMMFF94GPU();

      //! Assignment
      OBForceFieldMMFF94GPU &operator = (OBForceFieldMMFF94GPU &);

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


      // Helper function to get error string
      // *********************************************************************
      static const char* oclErrorString(int error)
      {
          static const char* errorString[] = {
              "CL_SUCCESS",
              "CL_DEVICE_NOT_FOUND",
              "CL_DEVICE_NOT_AVAILABLE",
              "CL_COMPILER_NOT_AVAILABLE",
              "CL_MEM_OBJECT_ALLOCATION_FAILURE",
              "CL_OUT_OF_RESOURCES",
              "CL_OUT_OF_HOST_MEMORY",
              "CL_PROFILING_INFO_NOT_AVAILABLE",
              "CL_MEM_COPY_OVERLAP",
              "CL_IMAGE_FORMAT_MISMATCH",
              "CL_IMAGE_FORMAT_NOT_SUPPORTED",
              "CL_BUILD_PROGRAM_FAILURE",
              "CL_MAP_FAILURE",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "",
              "CL_INVALID_VALUE",
              "CL_INVALID_DEVICE_TYPE",
              "CL_INVALID_PLATFORM",
              "CL_INVALID_DEVICE",
              "CL_INVALID_CONTEXT",
              "CL_INVALID_QUEUE_PROPERTIES",
              "CL_INVALID_COMMAND_QUEUE",
              "CL_INVALID_HOST_PTR",
              "CL_INVALID_MEM_OBJECT",
              "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
              "CL_INVALID_IMAGE_SIZE",
              "CL_INVALID_SAMPLER",
              "CL_INVALID_BINARY",
              "CL_INVALID_BUILD_OPTIONS",
              "CL_INVALID_PROGRAM",
              "CL_INVALID_PROGRAM_EXECUTABLE",
              "CL_INVALID_KERNEL_NAME",
              "CL_INVALID_KERNEL_DEFINITION",
              "CL_INVALID_KERNEL",
              "CL_INVALID_ARG_INDEX",
              "CL_INVALID_ARG_VALUE",
              "CL_INVALID_ARG_SIZE",
              "CL_INVALID_KERNEL_ARGS",
              "CL_INVALID_WORK_DIMENSION",
              "CL_INVALID_WORK_GROUP_SIZE",
              "CL_INVALID_WORK_ITEM_SIZE",
              "CL_INVALID_GLOBAL_OFFSET",
              "CL_INVALID_EVENT_WAIT_LIST",
              "CL_INVALID_EVENT",
              "CL_INVALID_OPERATION",
              "CL_INVALID_GL_OBJECT",
              "CL_INVALID_BUFFER_SIZE",
              "CL_INVALID_MIP_LEVEL",
              "CL_INVALID_GLOBAL_WORK_SIZE",
          };

          const int errorCount = sizeof(errorString) / sizeof(errorString[0]);

          const int index = -error;

          return (index >= 0 && index < errorCount) ? errorString[index] : "";

      }



  }; // class OBForceFieldMMFF94GPU

}// namespace OpenBabel

//! \file forcefieldmmff94gpu.h
//! \brief MMFF94 force field


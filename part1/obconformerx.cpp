/**********************************************************************
obconformer.cpp - Run a Monte Carlo conformer search for a molecule

Copyright (C) 2005-2007 Geoffrey R. Hutchison
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/obutil.h>

#include <openbabel/forcefield.h>

#include <stdio.h>
#include <iostream>
#include <fstream>

#ifdef _OPENMP // IF OPENMP SUPPORTED BY COMPILER
#include <omp.h> // in order to have access to OpenMP API (e.g. omp_get_thread_num())
#endif

#include <iostream>
#include <fstream>
#include "benchmark-utils.hpp"
#include <openbabel/mapkeys.h>


using namespace std;
using namespace OpenBabel;

OBRotorList rl;
OBRandom generator;

int main(int argc,char *argv[])
{
  if (argc != 4)
    {
      cout << "Usage: obconformer NSteps GeomSteps <file>" << endl;
      return(-1);
    }

  int weightSteps, geomSteps;
  weightSteps = atoi(argv[1]);
  geomSteps = atoi(argv[2]);

  string filename = argv[3];

  ifstream ifs(filename.c_str());
  if (!ifs)
    {
      cerr << "Error! Cannot read input file!" << endl;
      return(-1);
    }
  
  OBConversion conv(&ifs, &cout);
  OBFormat* pFormat;
  
  pFormat = conv.FormatFromExt(argv[3]);
  if ( pFormat == NULL )
    {
      cerr << "Error! Cannot read file format!" << endl;
      return(-1);
    }

  // Finally, we can do some work!
  OBMol mol;
  if (! conv.SetInAndOutFormats(pFormat, pFormat))
    {
      cerr << "Error! File format isn't loaded" << endl;
      return (-1);
    }

  string ff = "MMFF94"; // conformer search is only supported using MMFF94 forcefield

  OBForceField *pFF = OBForceField::FindForceField("MMFF94");
  pFF->SetLogFile(&cerr);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);

  // NEXT IS BOILER PLATE CODE USED TO INSTRUMENT THE APPLICATION IN ORDER TO
  // OBTAIN PERFORMANCE INFORMATION STATISTICS (RUNTIME AND MEMORY ALLOCATION)

  Timer totalTimer, readTimer, setupTimer, computeTimer;
  double totalTime, readTime, setupTime, computeTime;
  int nthreads = 1; // single core is just one thread
  unsigned int totalSteps = 1; // total steps is only used when doing minimization

	#ifdef _OPENMP // START OPENMP DEBUG
	#pragma omp parallel default(none) shared(nthreads)
	{
	  nthreads = omp_get_num_threads(); // determine the number of threads
	} // parallel region completes
	#endif // END OPENMP DEBUG

	// This code only deals with creating a file to store the timings for
	// the execution of the program along with other interesting statistics
	char cwd[1024]; // a reasonably large current working directory (cwd)
	char filepath[1100]; // a reasonably large file name plus cwd
	const char * statsfile = "obconformer_runstats";
	std::ofstream output;
	if (getcwd(cwd, sizeof (cwd)) != NULL) {
		std::cout << "Current working dir: " << cwd << std::endl;
	} else {
		std::cout << "The getcwd() method failed. Aborting." << std::endl;
		exit(1);
	}//if(getcwd()

	sprintf(filepath, "%s/%s_%s_t%d.mat", cwd, statsfile, ff.c_str(), nthreads);
	std::cout << "Writing output file to: " << filepath << std::endl;
	output.open(filepath, ios::out | ios::app ); // The file is open in append mode

	//              1     2       3        4        5          6          7        8         9          10         11    12
	output << "#METHOD: " << ff << " DATASET: " << filename.c_str() << std::endl;
	output << "#THREADS ENERGY MOL_MASS NUM_ATOMS NUM_ROTORS NUM_CONF TOT_TIME TIME_READ TIME_SETUP TIME_COMPUTE STEPS #MOL_NAME " << std::endl;

	// A second file is created to store information on the timings of specific parts
	// from the MMFF94 calculation routines. breakdown of time spent in calculations
	char filepath2[1100];
	std::ofstream output2;
	sprintf(filepath2, "%s/%s_%s_t%d_compute.mat", cwd, statsfile, ff.c_str(), nthreads);
	std::cout << "Writing method breakdown detail file to: " << filepath2 << std::endl;
	output2.open(filepath2, ios::out | ios::app ); // The file is open in append mode

  //           1      2       3        4         5     6     7      8       9         10         11      12
  output2 << "#METHOD: " << ff << " THREADS: " << nthreads << " DATASET: " << filename.c_str() << std::endl;
  output2 << "#E_BOND E_ANGLE E_STRBND E_TORSION E_OOP E_VDW E_ELEC N_ATOMS PAIRS_VDW PAIRS_ELEC MEM_VDW MEM_ELEC " << std::endl;

  double bondCalcTime, angleCalcTime, strbndCalcTime, torsionCalcTime, oopCalcTime, vdwCalcTime, electrostaticCalcTime;
  int numPairsVDW, numPairsElectrostatic;
    
  while(ifs.peek() != EOF && ifs.good()) {
    
    totalTimer.start();
    readTimer.start();
    
      mol.Clear();
      conv.Read(&mol);
      
      mol.AddHydrogens(); //NOTE: the original conformer search in OpenBabel doesn't have a switch for adding hydrogens

    readTime = readTimer.get();
    setupTimer.start();
      
      pFF->Setup(mol);
      
    setupTime = setupTimer.get();
    computeTimer.start();
      
      pFF->WeightedRotorSearch(weightSteps, geomSteps);
      pFF->ConjugateGradients(geomSteps); // final cleanup
      pFF->UpdateCoordinates(mol);
      
    computeTime = computeTimer.get();
    totalTime = totalTimer.get();
    
    // !! ENERGY CALCULATED HERE JUST TO COMPLETE LOGGING TABLE BUT NOT PART OF THE CONFORMER SEARCH IN OPEN BABEL
    double energy = pFF->Energy();

    // THREADS  ENERGY  MOL_MASS  NUM_ATOMS  NUM_ROTORS  NUM_CONF  TOT_TIME  TIME_READ  TIME_SETUP  TIME_COMPUTE STEPS  #MOL_NAME
    output << nthreads << " " << energy << " " << mol.GetExactMass() << " " << mol.NumAtoms()
    		<< " " << mol.NumRotors() << "  " <<  mol.NumConformers() << " "
    		<< totalTime <<  " " << readTime << " " << " " << setupTime << " " << computeTime << " "
    		<< totalSteps << " #" << mol.GetTitle()  // comment added to avoid errors when reading matrix in Octave
    		<< std::endl;

    map<string, double> timings = pFF->getTimings();
    map<string, size_t> memalloc = pFF->getAllocatedMemory();
    MapKeys mk;
    // 1      2       3        4         5     6     7      8       9         10         11      12
    // E_BOND E_ANGLE E_STRBND E_TORSION E_OOP E_VDW E_ELEC N_ATOMS PAIRS_VDW PAIRS_ELEC MEM_VDW MEM_ELEC
    output2 << timings[mk.TIME_BOND_CALCULATIONS] << " "  // 1
    		<< timings[mk.TIME_ANGLE_CALCULATIONS] << " " // 2
    		<< timings[mk.TIME_STRBND_CALCULATIONS] << " " // 3
    		<< timings[mk.TIME_TORSION_CALCULATIONS] << " " // 4
    		<< timings[mk.TIME_OOP_CALCULATIONS] << " " // 5
    		<< timings[mk.TIME_VDW_CALCULATIONS] << " " // 6
    		<< timings[mk.TIME_ELECTROSTATIC_CALCULATIONS] << " " // 7
    		<< mol.NumAtoms() << " " // 8
    		<< timings[mk.TOTAL_VDW_CALCULATIONS] << " " // 9
    		<< timings[mk.TOTAL_ELECTROSTATIC_CALCULATIONS] << " " // 10
    		<< memalloc[mk.MEM_VDW_CALCULATIONS] << " " // 11
    		<< memalloc[mk.MEM_ELECTROSTATIC_CALCULATIONS] << std::endl; // 12
      
      
      conv.Write(&mol);
    } // while reading molecules
    
    output.close();
    output2.close();
  
  return(0);
}

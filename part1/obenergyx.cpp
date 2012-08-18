/**********************************************************************
obenergy.cpp - calculate the energy for a molecule

Copyright (C) 2006 Tim Vandermeersch
 
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
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#ifndef _MSC_VER
  #include <unistd.h>
#endif

#ifdef _OPENMP // IF OPENMP SUPPORTED BY COMPILER
#include <omp.h> // in order to have access to OpenMP API (e.g. omp_get_thread_num())
#endif

#include <iostream>
#include <fstream>
#include "benchmark-utils.hpp"
#include <openbabel/mapkeys.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  int verbose = 0;
  bool hydrogens = false;
  string basename, filename = "", option, option2, ff = "";
  OBConversion conv;

  if (argc < 2) {
    cout << "Usage: obenergy [options] <filename>" << endl << endl;
    cout << "options:      description:" << endl << endl;
    cout << "  -v          verbose: print out indivual energy interactions" << endl << endl;
    cout << "  -h          add hydrogens before calculating energy" << endl << endl;
    cout << "  -ff ffid    select a forcefield" << endl << endl;
    cout << "              available forcefields:" << endl << endl;
    OBPlugin::List("forcefields", "verbose");
    exit(-1);
  } else {
    int ifile = 1;
    for (int i = 1; i < argc; i++) {
      option = argv[i];
      
      if (option == "-v") {
        verbose = 1;
        ifile++;
        break;
      }

      if (option == "-h") {
        hydrogens = true;
        ifile++;
      }

      if ((option == "-ff") && (argc > (i+1))) {
        ff = argv[i+1];
        ifile += 2;
      }
    }
    
    basename = filename = argv[ifile];
    size_t extPos = filename.rfind('.');
    size_t startPos = 0;

    if (extPos!= string::npos) {
      basename = filename.substr(startPos, extPos);
    }


  }


  // Find Input filetype
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());
    
  if (!format_in || !conv.SetInFormat(format_in)) {
    cerr << program_name << ": cannot read input format!" << endl;
    exit (-1);
  }

  ifstream ifs;
  ofstream ofs;

  // Read the file
  ifs.open(filename.c_str());
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }

  OBForceField* pFF = OBForceField::FindForceField(ff);
  if (!pFF) {
    cerr << program_name << ": could not find forcefield '" << ff << "'." <<endl;
    exit (-1);
  }
  pFF->SetLogFile(&cout);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);
//  if (verbose)
//    pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
//  else
//    pFF->SetLogLevel(OBFF_LOGLVL_MEDIUM);

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
	const char * statsfile = "obenergy_runstats";
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

	//              1     2       3        4        5          6          7        8         9          10        11    12
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

	// A third file is created to store information on the memory allocation of data structures
	// from the MMFF94 calculation routines. breakdown of memory allocated for each calculation type
	char filepath3[1100];
	std::ofstream output3;
	sprintf(filepath3, "%s/%s_%s_t%d_malloc.mat", cwd, statsfile, ff.c_str(), nthreads);
	std::cout << "Writing memory allocation breakdown detail file to: " << filepath3 << std::endl;
	output3.open(filepath3, ios::out | ios::app ); // The file is open in append mode

	//           1     2      3       4        5         6     7     8      9      10      11       12        13    14    15
	output3 << "#METHOD: " << ff << " THREADS: " << nthreads << " DATASET: " << filename.c_str() << std::endl;
	output3 << "#ATOMS M_BOND M_ANGLE M_STRBND M_TORSION M_OOP M_VDW M_ELEC C_BOND C_ANGLE C_STRBND C_TORSION C_OOP C_VDW C_ELEC" << std::endl;

  double bondCalcTime, angleCalcTime, strbndCalcTime, torsionCalcTime, oopCalcTime, vdwCalcTime, electrostaticCalcTime;
  int numPairsVDW, numPairsElectrostatic;

  OBMol mol;
  double energy;
  for (c=1;;c++) {
    mol.Clear();

    totalTimer.start();
    readTimer.start();

    if (!conv.Read(&mol, &ifs))
      break;
    if (mol.Empty())
      break;

    if (hydrogens)
      mol.AddHydrogens();
       
    readTime = readTimer.get();
    setupTimer.start();


    if (!pFF->Setup(mol)) {
      cerr << program_name << ": could not setup force field." << endl;
      exit (-1);
    }

    setupTime = setupTimer.get();
    computeTimer.start();
    
    energy = pFF->Energy(false);

    computeTime = computeTimer.get();
    totalTime = totalTimer.get();

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

	// 1     2      3       4        5         6     7     8      9      10      11       12        13    14    15
    // ATOMS M_BOND M_ANGLE M_STRBND M_TORSION M_OOP M_VDW M_ELEC C_BOND C_ANGLE C_STRBND C_TORSION C_OOP C_VDW C_ELEC
    output3 << mol.NumAtoms() << " " // 1
    		<< memalloc[mk.MEM_BOND_CALCULATIONS] << " " // 2
    		<< memalloc[mk.MEM_ANGLE_CALCULATIONS] << " " // 3
    		<< memalloc[mk.MEM_STRBND_CALCULATIONS] << " " // 4
    		<< memalloc[mk.MEM_TORSION_CALCULATIONS] << " " // 5
    		<< memalloc[mk.MEM_OOP_CALCULATIONS] << " " // 6
    		<< memalloc[mk.MEM_VDW_CALCULATIONS] << " " // 7
    		<< memalloc[mk.MEM_ELECTROSTATIC_CALCULATIONS] << " " // 8
    		<< timings[mk.TOTAL_BOND_CALCULATIONS] << " " // 9
    		<< timings[mk.TOTAL_ANGLE_CALCULATIONS] << " " // 10
    		<< timings[mk.TOTAL_STRBND_CALCULATIONS] << " " // 11
    		<< timings[mk.TOTAL_TORSION_CALCULATIONS] << " " // 12
    		<< timings[mk.TOTAL_OOP_CALCULATIONS] << " " // 13
    		<< timings[mk.TOTAL_VDW_CALCULATIONS] << " " // 14
    		<< timings[mk.TOTAL_ELECTROSTATIC_CALCULATIONS] << " " // 15
    		<< std::endl;

    if (!isfinite(energy)) {
      cerr << " Title: " << mol.GetTitle() << endl;
      FOR_ATOMS_OF_MOL(atom, mol) {
        cerr << " x: " << atom->x() << " y: " << atom->y() << " z: " << atom->z() << endl;
      }
    }

  } // end for loop

  output.close();
  output2.close();
  output3.close();

  return(1);
}

/* obenergy man page*/
/** \page calculate the energy for a molecule
*
* \n
* \par SYNOPSIS
*
* \b obenergy [options] \<filename\>
*
* \par DESCRIPTION
*
* The obenergy tool can be used to calculate the energy for molecules 
* inside (multi-)molecule files (e.g., MOL2, etc.)
*
* \par OPTIONS
*
* If no filename is given, obenergy will give all options including the
* available forcefields.
*
* \b -v:
*     Verbose: print out all individual energy interactions \n\n
* \b -ff \<forcefield\>:
*     Select the forcefield \n\n
*
* \par EXAMPLES
*  - View the possible options, including available forcefields: 
*   obenergy
*  - Calculate the energy for the molecule(s) in file test.mol2:
*   obenergy test.mol2
*  - Calculate the energy for the molecule(s) in file test.mol2 using the Ghemical forcefield:
*   obenergy -ff Ghemical test.mol2 
*  - Calculate the energy for the molecule(s) in file test.mol2 and print out all individual energy interactions:
*    obenergy -v test.mol2
*
* \par AUTHORS
*
* The obenergy program was contributed by \b Tim \b Vandermeersch.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 2007 by Tim Vandermeersch. \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at: http://openbabel.org/ \n
*   The web pages for Open Babel Molecular Mechanics can be found at: 
*   http://openbabel.org/wiki/Molecular_mechanics \n
**/

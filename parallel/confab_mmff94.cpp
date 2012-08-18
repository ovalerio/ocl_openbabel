/**********************************************************************
confab - Systematic conformer generation

Copyright (C) 2010 Noel O'Boyle
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#define CONFAB_VER "1.0.1"

#ifdef _OPENMP // START OPENMP DEBUG
#include <omp.h> // in order to have access to OpenMP API (e.g. omp_get_thread_num())
#endif

#include "hr_time.h" // this library is used to perform timings
#include <errno.h>

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/confsearch.h>

#ifdef _MSC_VER
	typedef char TCHAR;
	#include "getopt.h"
#else
	#include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;

class Confab
{
public:
  Confab(OBConversion &conv, double rmsd_cutoff, double energy_cutoff, int conf_cutoff, bool verbose, bool include_original);
  void DisplayConfig();
  void Run();
private:
  double rmsd_cutoff;
  double energy_cutoff;
  unsigned int conf_cutoff;
  bool verbose;
  bool include_original;
  OBConversion &conv;
  OBForceField *pff;
};

///////////////////////////////////////////////////////////////////////////////
//! \brief Find the molecule(s) with or without a given SMART pattern
int main(int argc,char **argv)
{
  char c;
  char *program_name = argv[0];
  char *iext;

  OBConversion conv;
  OBFormat *piFormat = NULL;
  OBFormat *poFormat = NULL;
  float rmsd_cutoff = 0.5, energy_cutoff = 50.0;
  unsigned int conf_cutoff = 1000000; // 1 Million
  bool verbose = false;
  bool include_original = false;

  // Parse options
  while ((c = getopt(argc, argv, "i:o:r:e:c:v:a")) != -1)
    {
#ifdef _WIN32
	    char optopt = c;
#endif
      switch (c)
        {
        case 'i':
          iext = optarg;

          piFormat = conv.FindFormat(iext);
	  
          if(piFormat==NULL)
            {
              cerr << program_name << ": cannot read input format!" << endl;
              exit(-1);
            }
          break;
        case 'o':
          iext = optarg;

          poFormat = conv.FindFormat(iext);
	  
          if(poFormat==NULL)
            {
              cerr << program_name << ": cannot read output format!" << endl;
              exit(-1);
            }
          break;
        case 'r': // set rmsd_cutoff

          c = sscanf(optarg, "%f", &rmsd_cutoff);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -r option" << endl;
              exit (-1);
            }
          break;
        case 'e': // set energy_cutoff

          c = sscanf(optarg, "%f", &energy_cutoff);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -e option" << endl;
              exit (-1);
            }
          break;
        case 'c': // set max number of conformations

          c = sscanf(optarg, "%i", &conf_cutoff);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -c option" << endl;
              exit (-1);
            }
          break;
        case 'v': // verbose
          verbose = true;
          break;
        case 'a': // include input structure as first conformer
          include_original = true;
          break;
        case '?':
          if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf (stderr,
                     "Unknown option character `\\x%x'.\n",
                     optopt);
          return 1;
        }
    }
  int index = optind;

  if (argc-index != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " [options] <inputfile> <outputfile>\n";
      err += "Options:\n";
      err += "   -i <format> Input file format\n";
      err += "   -o <format> Output file format\n";
      err += "   -r <rmsd>   RMSD cutoff (default 0.5 Angstrom)\n";
      err += "   -e <energy> Energy cutoff (default 50.0 kcal/mol)\n";
      err += "   -c <#confs> Max number of conformers to test (default is 1 million)\n";
      err += "   -a          Include the input conformation as the first conformer\n";
      err += "   -v          Verbose - display information on torsions found\n";
      cerr << err << ends;
      exit(-1);
    }

  // Set up the input and output files
  char* inputfile = argv[index++];
  char* outputfile = argv[index];

  ifstream ifs;
  ifs.open(inputfile);
  if (!ifs) {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
  }

  // Find Input/Output format (if not already set)
  if (!piFormat) {
    piFormat = conv.FormatFromExt(inputfile);
    if (!piFormat) {
      cerr << program_name << ": cannot read input format!" << endl;
      return (-1);
    }
  }
  if (!poFormat) {
    poFormat = conv.FormatFromExt(outputfile);
    if (!poFormat) {
      cerr << program_name << ": cannot read output format!" << endl;
      return (-1);
    }
  }

  ofstream ofs;
  ofs.open(outputfile);
  if (!ifs) {
      cerr << program_name << ": cannot write output file!" << endl;
      exit (-1);
  }

  // Run Confab!!
  conv.SetInStream(&ifs);
  conv.SetOutStream(&ofs);
  conv.SetInAndOutFormats(piFormat, poFormat);

  cout << "**Starting Confab " << CONFAB_VER << endl;
  cout << "..Input file = " << inputfile << endl;
  cout << "..Output file = " << outputfile << endl;

  Confab confab(conv, rmsd_cutoff, energy_cutoff, conf_cutoff, verbose, include_original);
  confab.DisplayConfig();
  confab.Run();

  cout << endl << "**Finishing Confab" << endl;

}

Confab::Confab(OBConversion &conv, double rmsd_cutoff, double energy_cutoff, int conf_cutoff, bool verbose, bool include_original):
conv(conv), rmsd_cutoff(rmsd_cutoff), energy_cutoff(energy_cutoff), conf_cutoff(conf_cutoff), verbose(verbose), include_original(include_original)
{
  pff = OpenBabel::OBForceField::FindType("mmff94");
  if (!pff) {
    cout << "Cannot find forcefield!" << endl;
    exit(-1);
  }
  if (verbose) {
    pff->SetLogLevel(OBFF_LOGLVL_LOWER);
    pff->SetLogFile(&cout);
  }
}

void Confab::Run()
{

	OBMol mol;
	unsigned int N = 0;
	bool moleculeRead = false;

	// some variables needed to compute statistics on the code
	int nthreads = 1; // the number of processors
	CStopWatch stopWatchTotal; // this class will be used in case no OpenMP timer available
	CStopWatch stopWatchCompute; // only the computationally intensive part is measured here
	double startTime, endTime, elapsedTime; // elapsedTime execution time in seconds
	double startCompute, endCompute, elapsedCompute; // consider only the computation time

#ifdef _OPENMP // START OPENMP DEBUG
#pragma omp parallel default(none) shared(nthreads)
	{
		nthreads = omp_get_num_threads();
	}
#endif

	// This code only deals with creating a file to store the timings for
	// the execution of the program along with other interesting statistics
	char cwd[1024]; // a reasonably large current working directory (cwd)
	char filepath[1100]; // a reasonably large file name plus cwd
	char * filename = "confab_runstats";
	FILE * output;
	if (getcwd(cwd, sizeof (cwd)) != NULL) {
		fprintf(stdout, "Current working dir: %s\n", cwd);
	} else {
		perror("The getcwd() method failed. Aborting.\n");
		exit(1);
	}//if(getcwd()

	sprintf(filepath, "%s/%s.dat", cwd, filename);
	printf("Writing output file to: %s\n", filepath);
	output = fopen(filepath, "a"); // The file is open in append mode

	//	fprintf(output, "NTHREADS MOLNUM TITLE BONDS CONFORMERS TOTTIME CPUTIME CPU/TOT CONF/CPUTIME \n");
	fprintf(output, "#P   MOL TITLE    ROT  CONF        TIME      CPU_TIME    CPU/TOT  CONF/CPUTIME \n");

	while (true) { // loop processing molecules from the file until the file is finished

#ifdef _OPENMP // START OPENMP STOPWATCH
		startTime = omp_get_wtime();
		startCompute = omp_get_wtime();
#else // START UNIX STOPWATCH
		stopWatchTotal.startTimer();
		stopWatchCompute.startTimer();
#endif

		moleculeRead = conv.Read(&mol); // find if there is another molecule in the file
		if (!moleculeRead) break; // break out from the loop
		N += 1; // increase the molecule counter N

		cout << "**Molecule " << N << endl << "..title = " << mol.GetTitle() << endl;
		cout << "..number of rotatable bonds = " << mol.NumRotors() << endl;
		mol.AddHydrogens();
		bool success = pff->Setup(mol);
		if (!success) {
			cout << "!!Cannot set up forcefield for this molecule\n"
					<< "!!Skipping\n" << endl;
			continue;
		}
		pff->DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff);

		pff->GetConformers(mol);

#ifdef _OPENMP // START OPENMP STOPWATCH
		endCompute = omp_get_wtime();
		elapsedCompute = endCompute - startCompute;
#else // START UNIX STOPWATCH
		stopWatchCompute.stopTimer();
		elapsedCompute = stopWatchCompute.getElapsedTime();
#endif

		int nconfs = include_original ? mol.NumConformers() : mol.NumConformers() - 1;

		cout << "..generated " << nconfs << " conformers" << endl;

		unsigned int c = include_original ? 0 : 1;

		for (; c < mol.NumConformers(); ++c) {
			mol.SetConformer(c);
			if(!conv.GetOutFormat()->WriteMolecule(&mol, &conv))
				break;
		}
		cout << endl;

#ifdef _OPENMP // STOP OPENMP STOPWATCH
		endTime = omp_get_wtime();
		elapsedTime = endTime - startTime;
#else // STOP UNIX STOPWATCH
		stopWatchTotal.stopTimer();
		elapsedTime = stopWatchTotal.getElapsedTime();
#endif

		printf("Threads = %d  Elapsed time = %f\n", nthreads, elapsedTime);

		fprintf(output, " %-3d %-3d %.8s %-3d  %-9d  %-9.4f  %-10.4f  %-10.4f  %-10.4f \n",
				nthreads, N, mol.GetTitle(), mol.NumRotors(), mol.NumConformers(),
				elapsedTime, elapsedCompute, elapsedCompute/elapsedTime,
				mol.NumConformers()/elapsedCompute);
		fflush(output);

	}//while(moreMolecules)

	fclose(output);

}//Run()

void Confab::DisplayConfig()
{
  cout << "..Input format = " << conv.GetInFormat()->GetID() << endl;
  cout << "..Output format = " << conv.GetOutFormat()->GetID() << endl;
  cout << "..RMSD cutoff = " << rmsd_cutoff << endl;
  cout << "..Energy cutoff = " << energy_cutoff << endl;
  cout << "..Conformer cutoff = " << conf_cutoff << endl;
  cout << "..Write input conformation? " << (include_original ? "True" : "False") << endl;
  cout << endl;
}

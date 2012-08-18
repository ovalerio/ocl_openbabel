/*
 * A minimal validation test for the OpenBabel MMFF94 enery computation.
 * It makes use of the forcefields available from OpenBabel project to execute.
 */
#include <openbabel/babelconfig.h>

#include <fstream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

void ComputeEnergies(string molecules_file, string results_file)
{
	std::ifstream ifs;
	if (!SafeOpen(ifs, molecules_file.c_str()))
		return;

	std::ofstream ofs;
	if (!SafeOpen(ofs, results_file.c_str()))
		return;

	OBMol mol;
	OBConversion conv(&ifs, &cout);
	char buffer[BUFF_SIZE];

	if(! conv.SetInAndOutFormats("SDF","SDF")){
		cerr << "SDF format is not loaded" << endl;
		return;
	}

	OBForceField* pFF = OBForceField::FindForceField("MMFF94");

	if (pFF == NULL) {
		cerr << "Cannot load force field!" << endl;
		return;
	}

	pFF->SetLogFile(&cout);
	pFF->SetLogLevel(OBFF_LOGLVL_NONE);

	double energy;
	while (ifs)
	{
		mol.Clear();
		conv.Read(&mol);
		if (mol.Empty())
			continue;

		if (!pFF->Setup(mol)) {
			cerr << "Could not setup force field on molecule: " << mol.GetTitle() << endl;
			return;
		}

		// Don't compute gradients
		// Write the energy computation to the file
//		energy = pFF->Energy(false);
		energy = pFF->Energy(false);
		sprintf(buffer, "Molecule: %s - Calc. Energy: %15.5f\n", mol.GetTitle(), pFF->Energy(false));
		ofs << buffer;
	}

	cout << " MMFF94 force field energies written successfully" << endl;
	return;
}

int main(int argc,char *argv[])
{
	// turn off slow sync with C-style output (we don't use it anyway).
	std::ios::sync_with_stdio(false);

	// Define location of file formats for testing
#ifdef FORMATDIR
	char env[BUFF_SIZE];
	snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
	putenv(env);
#endif

	string molecules_file = "forcefield.sdf";
	string results_file = "results.txt";

	if (argc > 3)
	{
		molecules_file = argv[2];
		results_file = argv[3];
	}

	ComputeEnergies(molecules_file, results_file);
	return 0;

}

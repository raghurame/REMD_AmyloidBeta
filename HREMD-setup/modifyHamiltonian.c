#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/*
Segments in topol.top

[ moleculetype ]
[ atoms ]
[ bonds ]
[ pairs ]
[ angles ]
[ dihedrals ]
[ position_restraints ]
[ system ]
[ molecules ]
*/

typedef struct topologyBool
{
	int moleculeType, atoms, bonds, pairs, angles, dihedrals, posre, system, molecules;
} TOPOLOGY_BOOL;

typedef struct topologyAtoms
{
	int nr, resnr, cgnr;
	char type[10], residue[10], atom[10];
	float charge, mass;
} TOPOLOGY_ATOMS;

typedef struct topologyBonds
{
	int ai, aj, func;
	char includeString[10];
} TOPOLOGY_BONDS;

typedef struct topologyPairs
{
	int ai, aj, func;
} TOPOLOGY_PAIRS;

typedef struct topologyAngles
{
	int ai, aj, ak, func;
	char includeString[10];
} TOPOLOGY_ANGLES;

typedef struct topologyDihedrals
{
	int ai, aj, ak, al, func;
	char includeString[10];
} TOPOLOGY_DIHEDRALS;

typedef struct posre
{
	int i, func, fcx, fcy, fcz;
} POSITIONAL_RESTRAINTS;

// System directive can occur multiple times in the top file
typedef struct system
{
	char name[1000];
} TOPOLOGY_SYSTEM;

// Molecule directive can occur multiple times in the top file
typedef struct molecule
{
	char name[1000];
	int nMolecules;
} TOPOLOGY_MOLECULE;

typedef struct bondedDefines
{
	char defineString, defineType;
	float value, constant;
	int np;
} BONDED_DEFINES;

TOPOLOGY_ATOMS *readTopAtoms (FILE *topolTopITP, TOPOLOGY_BOOL topCurrentPosition, TOPOLOGY_ATOMS *inputAtoms, int *nAtoms)
{
	rewind (topolTopITP);
	char lineString[2000], atomString[2000];

	(*nAtoms) = 0;

	// Counting the number of entries in the [ atoms ] directive
	while (fgets (lineString, 2000, topolTopITP) != NULL)
	{
		if (lineString[0] != ';')
		{
			if ((topCurrentPosition.atoms == 1) && (lineString[0] == '[')) {
				topCurrentPosition.atoms = 0; }

			if (topCurrentPosition.atoms == 1)
			{
				for (int i = 0; lineString[i] != ';'; ++i) {
					atomString[i] = lineString[i];
					atomString[i + 1] = '\0'; }

				(*nAtoms)++;
			}

			if (strstr (lineString, "[ atoms ]")) {
				topCurrentPosition.atoms = 1; }
		}
	}

	// One empty line was also counted in the previous loop
	printf("Number of atoms detected in topology file: %d\n", (*nAtoms) - 1);

	inputAtoms = (TOPOLOGY_ATOMS *) malloc ((*nAtoms) * sizeof (TOPOLOGY_ATOMS));

	rewind (topolTopITP);
	int currentAtom = 0;
	// Storing the information under [ atoms ] directive
	while (fgets (lineString, 2000, topolTopITP) != NULL)
	{
		if (lineString[0] != ';')
		{
			if ((topCurrentPosition.atoms == 1) && (lineString[0] == '[')) {
				topCurrentPosition.atoms = 0; }

			if (topCurrentPosition.atoms == 1)
			{
				for (int i = 0; lineString[i] != ';'; ++i) {
					atomString[i] = lineString[i];
					atomString[i + 1] = '\0'; }

				sscanf (atomString, "%d\n", &currentAtom);
				sscanf (atomString, "%d %s %d %s %s %d %f %f", &inputAtoms[currentAtom - 1].nr, &inputAtoms[currentAtom - 1].type, &inputAtoms[currentAtom - 1].resnr, &inputAtoms[currentAtom - 1].residue, &inputAtoms[currentAtom - 1].atom, &inputAtoms[currentAtom - 1].cgnr, &inputAtoms[currentAtom - 1].charge, &inputAtoms[currentAtom - 1].mass);
			}

			if (strstr (lineString, "[ atoms ]")) {
				topCurrentPosition.atoms = 1; }
		}
	}

	// To compensate for the extra empty line counted in the first 'while' loop
	(*nAtoms) -= 1;

	return inputAtoms;
}

BONDED_DEFINES *readBondedITP (FILE *ffBondedITP, TOPOLOGY_BOOL topCurrentPosition, BONDED_DEFINES **inputBondDefines, BONDED_DEFINES **inputAngleDefines, BONDED_DEFINES **inputProperDihedralDefines, BONDED_DEFINES **inputImproperDihedralDefines, int *nBondDefines, int *nAngleDefines, int *nProperDihedralDefines, int *nImproperDihedralDefines)
{
	rewind (ffBondedITP);
	char lineString[2000];

	// Counting the number of define statements corresponding to bonds, angles, dihedrals, and impropers
	(*nBondDefines) = 0;
	(*nAngleDefines) = 0;
	(*nProperDihedralDefines) = 0;
	(*nImproperDihedralDefines) = 0;

	while (fgets (lineString, 2000, ffBondedITP) != NULL)
	{
		if (lineString[0] != ';' && lineString[0] == '#')
		{
			if (strstr (lineString, "gb_")) {
				(*nBondDefines)++; }

			if (strstr (lineString, "ga_")) {
				(*nAngleDefines)++; }

			if (strstr (lineString, "gi_")) {
				(*nImproperDihedralDefines)++; }

			if (strstr (lineString, "gd_")) {
				(*nProperDihedralDefines)++; }
		}
	}

	printf("N defines:\n\n  Bonds: %d\n  Angles: %d\n  Proper dihedrals: %d\n  Improper dihedrals: %d\n\n", (*nBondDefines), (*nAngleDefines), (*nProperDihedralDefines), (*nImproperDihedralDefines));
}

int main(int argc, char const *argv[])
{
	FILE *ffBondedITP, *ffNonbondedITP, *topolTopITP, *ffBondedITP_output, *ffNonbondedITP_output, *topolTopITP_output;
	ffBondedITP = fopen (argv[1], "r");
	ffNonbondedITP = fopen (argv[2], "r");
	topolTopITP = fopen (argv[3], "r");
	ffBondedITP_output = fopen (argv[4], "w");
	ffNonbondedITP_output = fopen (argv[5], "w");
	topolTopITP_output = fopen (argv[6], "w");

	/*
		Ignore the lines starting with ';'
		Always take the lines starting with '#'
	*/

	TOPOLOGY_BOOL topCurrentPosition;
	TOPOLOGY_ATOMS *inputAtoms;
	TOPOLOGY_BONDS *inputBonds;
	TOPOLOGY_PAIRS *inputPairs;
	TOPOLOGY_ANGLES *inputAngles;
	TOPOLOGY_DIHEDRALS *inputDihedrals;
	POSITIONAL_RESTRAINTS *inputRestraints;
	TOPOLOGY_SYSTEM *inputSystem;
	TOPOLOGY_MOLECULE *inputMolecule;
	BONDED_DEFINES *inputBondDefines, *inputAngleDefines, *inputProperDihedralDefines, *inputImproperDihedralDefines;

	int nAtoms, nBondDefines, nAngleDefines, nProperDihedralDefines, nImproperDihedralDefines;

	inputAtoms = readTopAtoms (topolTopITP, topCurrentPosition, inputAtoms, &nAtoms);
	readBondedITP (ffBondedITP, topCurrentPosition, &inputBondDefines, &inputAngleDefines, &inputProperDihedralDefines, &inputImproperDihedralDefines, &nBondDefines, &nAngleDefines, &nProperDihedralDefines, &nImproperDihedralDefines);

	fclose (ffBondedITP);
	fclose (ffNonbondedITP);
	fclose (topolTopITP);
	fclose (ffBondedITP_output);
	fclose (ffNonbondedITP_output);
	fclose (topolTopITP_output);
	return 0;
}
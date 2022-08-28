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
} SYSTEM;

// Molecule directive can occur multiple times in the top file
typedef struct molecule
{
	char name[1000];
	int nMolecules;
} MOLECULE;

TOPOLOGY_ATOMS *readTopAtoms (FILE *topolTopITP, TOPOLOGY_BOOL topCurrentPosition, TOPOLOGY_ATOMS *inputAtoms)
{
	rewind (topolTopITP);
	char lineString[2000];

	while (fgets (lineString, 2000, topolTopITP) != NULL)
	{
			// Count the number of directives, include statements in this manner. If necessary, make a new function for this purpose.
		if (lineString[0] != ';')
		{
			if (lineString[0] == '#' || lineString[0] == '[')
			{
				printf("%s", lineString);
				sleep (1);
			}
		}
	}
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
	SYSTEM *inputSystem;
	MOLECULE *inputMolecule;

	inputAtoms = readTopAtoms (topolTopITP, topCurrentPosition, inputAtoms);

	fclose (ffBondedITP);
	fclose (ffNonbondedITP);
	fclose (topolTopITP);
	fclose (ffBondedITP_output);
	fclose (ffNonbondedITP_output);
	fclose (topolTopITP_output);
	return 0;
}
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
	int moleculeType, atoms, bonds, pairs, angles, dihedrals, posre, system, molecules, atomTypes, nonbondedParams;
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
	char defineString[10], defineType[10];
	float value, constant;
	int np;
} BONDED_DEFINES;

typedef struct nonbondedAtomtypes
{
	char name[10], ptype[10];
	int atomicNumber;
	float atomicMass, atomicCharge, c6, c12;
} NONBONDED_ATOMTYPES;

typedef struct nonbondedParams
{
	char i[10], j[10];
	int func;
	float c6, c12;
} NONBONDED_PARAMS;

typedef struct hotResidues
{
	char resName[10];
} HOT_RESIDUES;

typedef struct hotInteractions
{
	char resName1[10], resName2[10];
} HOT_INTERACTIONS;

TOPOLOGY_ATOMS *readTopAtoms (FILE *topolTopITP, TOPOLOGY_BOOL topCurrentPosition, TOPOLOGY_ATOMS *inputAtoms, int *nAtoms, FILE *topolTopITP_output, float lambda)
{
	rewind (topolTopITP);
	char lineString[2000], atomString[2000];

	// Resetting the position bool
	topCurrentPosition.atoms = 0;

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

	// Resetting the position bool
	topCurrentPosition.atoms = 0;

	// To compensate for the extra empty line counted in the first 'while' loop
	(*nAtoms) -= 1;

	// One empty line was also counted in the previous loop
	printf("Number of atoms detected in topology file: %d\n", (*nAtoms));

	inputAtoms = (TOPOLOGY_ATOMS *) malloc ((*nAtoms + 1) * sizeof (TOPOLOGY_ATOMS));

	rewind (topolTopITP);
	int currentAtom = 0;

	// Storing the information under [ atoms ] directive
	while (fgets (lineString, 2000, topolTopITP) != NULL)
	{
		if (lineString[0] != ';')
		{
			if ((topCurrentPosition.atoms == 1) && (lineString[0] == '[')) {
				topCurrentPosition.atoms = 0;
				fprintf(topolTopITP_output, "\n"); }

			if (topCurrentPosition.atoms == 1)
			{
				for (int i = 0; lineString[i] != ';'; ++i) {
					atomString[i] = lineString[i];
					atomString[i + 1] = '\0'; }

				sscanf (atomString, "%d\n", &currentAtom);
				sscanf (atomString, "%d %s %d %s %s %d %f %f", &inputAtoms[currentAtom - 1].nr, &inputAtoms[currentAtom - 1].type, &inputAtoms[currentAtom - 1].resnr, &inputAtoms[currentAtom - 1].residue, &inputAtoms[currentAtom - 1].atom, &inputAtoms[currentAtom - 1].cgnr, &inputAtoms[currentAtom - 1].charge, &inputAtoms[currentAtom - 1].mass);

				fprintf (topolTopITP_output, "%d\t\t%s\t\t%d\t\t%s\t\t%s\t\t%d\t\t%f\t\t%f\n", inputAtoms[currentAtom - 1].nr, inputAtoms[currentAtom - 1].type, inputAtoms[currentAtom - 1].resnr, inputAtoms[currentAtom - 1].residue, inputAtoms[currentAtom - 1].atom, inputAtoms[currentAtom - 1].cgnr, inputAtoms[currentAtom - 1].charge * lambda, inputAtoms[currentAtom - 1].mass);
			}

			if (strstr (lineString, "[ atoms ]")) {
				fprintf(topolTopITP_output, "%s", lineString);
				topCurrentPosition.atoms = 1; }

			if (topCurrentPosition.atoms == 0) {
				fprintf(topolTopITP_output, "%s", lineString); }
		}
	}

	return inputAtoms;
}

void readBondedITP (FILE *ffBondedITP, TOPOLOGY_BOOL topCurrentPosition, BONDED_DEFINES **inputBondDefines, BONDED_DEFINES **inputAngleDefines, BONDED_DEFINES **inputProperDihedralDefines, BONDED_DEFINES **inputImproperDihedralDefines, int *nBondDefines, int *nAngleDefines, int *nProperDihedralDefines, int *nImproperDihedralDefines, FILE *ffBondedITP_output, float lambda)
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

	// Store the values in define statements
	(*inputBondDefines) = (BONDED_DEFINES *) malloc ((*nBondDefines) * sizeof (BONDED_DEFINES));
	(*inputAngleDefines) = (BONDED_DEFINES *) malloc ((*nAngleDefines) * sizeof (BONDED_DEFINES));
	(*inputProperDihedralDefines) = (BONDED_DEFINES *) malloc ((*nProperDihedralDefines) * sizeof (BONDED_DEFINES));
	(*inputImproperDihedralDefines) = (BONDED_DEFINES *) malloc ((*nImproperDihedralDefines) * sizeof (BONDED_DEFINES));

	rewind (ffBondedITP);

	int currentBondDefine = 0, currentAngleDefine = 0, currentProperDihedralDefine = 0, currentImproperDihedralDefine = 0;
	char *defineType_local;
	defineType_local = (char *) malloc (500 * sizeof (char));

	while (fgets (lineString, 2000, ffBondedITP) != NULL)
	{
		if (lineString[0] != ';' && lineString[0] == '#')
		{
			if (strstr (lineString, "gb_")) {
				sscanf (lineString, "%s %s %f %f\n", 
					&(*inputBondDefines)[currentBondDefine].defineString,
					&(*inputBondDefines)[currentBondDefine].defineType, 
					&(*inputBondDefines)[currentBondDefine].value, 
					&(*inputBondDefines)[currentBondDefine].constant);

				fprintf(ffBondedITP_output, "%s\t%s\t%f\t%f\n", 
					(*inputBondDefines)[currentBondDefine].defineString,
					(*inputBondDefines)[currentBondDefine].defineType, 
					(*inputBondDefines)[currentBondDefine].value, 
					(*inputBondDefines)[currentBondDefine].constant * lambda);

				currentBondDefine++; }

			if (strstr (lineString, "ga_")) {
				sscanf (lineString, "%s %s %f %f\n", 
					&(*inputAngleDefines)[currentAngleDefine].defineString, 
					&(*inputAngleDefines)[currentAngleDefine].defineType, 
					&(*inputAngleDefines)[currentAngleDefine].value, 
					&(*inputAngleDefines)[currentAngleDefine].constant);

				fprintf(ffBondedITP_output, "%s\t%s\t%f\t%f\n", 
					(*inputAngleDefines)[currentAngleDefine].defineString, 
					(*inputAngleDefines)[currentAngleDefine].defineType, 
					(*inputAngleDefines)[currentAngleDefine].value, 
					(*inputAngleDefines)[currentAngleDefine].constant * lambda);

				currentAngleDefine++; }

			if (strstr (lineString, "gi_")) {
				sscanf (lineString, "%s %s %f %f\n", 
					&(*inputImproperDihedralDefines)[currentImproperDihedralDefine].defineString, 
					&(*inputImproperDihedralDefines)[currentImproperDihedralDefine].defineType, 
					&(*inputImproperDihedralDefines)[currentImproperDihedralDefine].value, 
					&(*inputImproperDihedralDefines)[currentImproperDihedralDefine].constant);

				fprintf(ffBondedITP_output, "%s\t%s\t%f\t%f\n", 
					(*inputImproperDihedralDefines)[currentImproperDihedralDefine].defineString, 
					(*inputImproperDihedralDefines)[currentImproperDihedralDefine].defineType, 
					(*inputImproperDihedralDefines)[currentImproperDihedralDefine].value, 
					(*inputImproperDihedralDefines)[currentImproperDihedralDefine].constant * lambda);

				currentImproperDihedralDefine++; }

			if (strstr (lineString, "gd_")) {
				sscanf (lineString, "%s %s %f %f %d\n", 
					&(*inputProperDihedralDefines)[currentProperDihedralDefine].defineString, 
					&(*inputProperDihedralDefines)[currentProperDihedralDefine].defineType, 
					&(*inputProperDihedralDefines)[currentProperDihedralDefine].value, 
					&(*inputProperDihedralDefines)[currentProperDihedralDefine].constant, 
					&(*inputProperDihedralDefines)[currentProperDihedralDefine].np);

				fprintf(ffBondedITP_output, "%s\t%s\t%f\t%f\t%d\n", 
					(*inputProperDihedralDefines)[currentProperDihedralDefine].defineString, 
					(*inputProperDihedralDefines)[currentProperDihedralDefine].defineType, 
					(*inputProperDihedralDefines)[currentProperDihedralDefine].value, 
					(*inputProperDihedralDefines)[currentProperDihedralDefine].constant * lambda, 
					(*inputProperDihedralDefines)[currentProperDihedralDefine].np);

				currentProperDihedralDefine++; }
		}
	}
}

void readNonbondedITP (FILE *ffNonbondedITP, TOPOLOGY_BOOL topCurrentPosition, NONBONDED_ATOMTYPES **inputNonbondedAtomtypes, NONBONDED_PARAMS **inputNonbondedParams, int *nNonbondedAtomtypes, int *nNonbondedParams, FILE *ffNonbondedITP_output, float lambda)
{
	char lineString[2000], rawString[2000];

	redocounting: ;
	topCurrentPosition.atomTypes = 0;
	topCurrentPosition.nonbondedParams = 0;
	(*nNonbondedAtomtypes) = 0;
	(*nNonbondedParams) = 0;

	int nNonbondedAtomtypes_local = 0, nNonbondedParams_local = 0;

	rewind (ffNonbondedITP);

	// Count the number of entries under each directive
	while (fgets (rawString, 2000, ffNonbondedITP) != NULL)
	{
		if (rawString[0] != ';')
		{
			for (int i = 0; rawString[i] != ';'; ++i) {
				lineString[i] = rawString[i];
				lineString[i + 1] = '\0'; }

			if (strlen (lineString) > 1)
			{				
				if (topCurrentPosition.atomTypes == 1 && lineString[0] == '[') 	{
					topCurrentPosition.atomTypes = 0; }

				if (topCurrentPosition.atomTypes == 1) {
					nNonbondedAtomtypes_local++; }

				if (strstr (lineString, "[ atomtypes ]")) {
					topCurrentPosition.atomTypes = 1; }

				if (topCurrentPosition.nonbondedParams == 1 && lineString[0] == '[') {
					topCurrentPosition.nonbondedParams = 0; }

				if (topCurrentPosition.nonbondedParams == 1) {
					nNonbondedParams_local++; }

				if (strstr (lineString, "[ nonbond_params ]")) {
					topCurrentPosition.nonbondedParams = 1; }
			}
		}
	}

	if (nNonbondedAtomtypes_local == 0 || nNonbondedParams_local == 0) {
		printf("Counting again...\n");
		fflush (stdout);
		sleep (10);
		goto redocounting; }

	printf("\nAllocating %d memory for (*inputNonbondedAtomtypes)\nAllocating %d memory for (*inputNonbondedParams)\n\n", nNonbondedAtomtypes_local, nNonbondedParams_local);

	// Allocating memory
	(*inputNonbondedAtomtypes) = (NONBONDED_ATOMTYPES *) malloc (nNonbondedAtomtypes_local * sizeof (NONBONDED_ATOMTYPES));
	(*inputNonbondedParams) = (NONBONDED_PARAMS *) malloc (nNonbondedParams_local * sizeof (NONBONDED_PARAMS));

	(*nNonbondedAtomtypes) = nNonbondedAtomtypes_local;
	(*nNonbondedParams) = nNonbondedParams_local;

	topCurrentPosition.atomTypes = 0;
	topCurrentPosition.nonbondedParams = 0;

	int currentAtomtype = 0, currentParams = 0;

	// Save information
	rewind (ffNonbondedITP);
	while (fgets (rawString, 2000, ffNonbondedITP) != NULL)
	{
		if (rawString[0] != ';')
		{
			for (int i = 0; rawString[i] != ';'; ++i) {
				lineString[i] = rawString[i];
				lineString[i + 1] = '\0'; }

			if (topCurrentPosition.atomTypes == 1 && lineString[0] == '[') 	{
				topCurrentPosition.atomTypes = 0;
				fprintf(ffNonbondedITP_output, "%s\n", lineString); }

			if (topCurrentPosition.atomTypes == 1)
			{
				sscanf (lineString, "%s %d %f %f %s %f %f\n", 
					&(*inputNonbondedAtomtypes)[currentAtomtype].name, 
					&(*inputNonbondedAtomtypes)[currentAtomtype].atomicNumber, 
					&(*inputNonbondedAtomtypes)[currentAtomtype].atomicMass, 
					&(*inputNonbondedAtomtypes)[currentAtomtype].atomicCharge, 
					&(*inputNonbondedAtomtypes)[currentAtomtype].ptype, 
					&(*inputNonbondedAtomtypes)[currentAtomtype].c6, 
					&(*inputNonbondedAtomtypes)[currentAtomtype].c12);

				fprintf(ffNonbondedITP_output, "%s %d %f %f %s %f %f\n", 
					(*inputNonbondedAtomtypes)[currentAtomtype].name, 
					(*inputNonbondedAtomtypes)[currentAtomtype].atomicNumber, 
					(*inputNonbondedAtomtypes)[currentAtomtype].atomicMass, 
					(*inputNonbondedAtomtypes)[currentAtomtype].atomicCharge, 
					(*inputNonbondedAtomtypes)[currentAtomtype].ptype, 
					(*inputNonbondedAtomtypes)[currentAtomtype].c6 * lambda, 
					(*inputNonbondedAtomtypes)[currentAtomtype].c12 * lambda);

				currentAtomtype++;
			}

			if (strstr (lineString, "[ atomtypes ]")) {
				topCurrentPosition.atomTypes = 1;
				fprintf(ffNonbondedITP_output, "%s\n", lineString); }

			if (topCurrentPosition.nonbondedParams == 1 && lineString[0] == '[') {
				topCurrentPosition.nonbondedParams = 0;
				fprintf(ffNonbondedITP_output, "%s\n", lineString); }

			if (topCurrentPosition.nonbondedParams == 1)
			{
				sscanf (lineString, "%s %s %d %f %f\n", 
					&(*inputNonbondedParams)[currentParams].i, 
					&(*inputNonbondedParams)[currentParams].j, 
					&(*inputNonbondedParams)[currentParams].func, 
					&(*inputNonbondedParams)[currentParams].c6, 
					&(*inputNonbondedParams)[currentParams].c12);

				currentParams++;
			}

			if (strstr (lineString, "[ nonbond_params ]")) {
				topCurrentPosition.nonbondedParams = 1; }
		}
	}
}

HOT_RESIDUES *readHotResNames (FILE *hotResidues, HOT_RESIDUES *hotResNames, int *nHotResNames)
{
	// Reading the number of lines in input file
	int nLines = 0;
	char lineString[1000];

	while (fgets (lineString, 1000, hotResidues) != NULL) {
		nLines++; }

	hotResNames = (HOT_RESIDUES *) malloc (nLines * sizeof (HOT_RESIDUES));
	(*nHotResNames) = nLines;

	// Storing the values from input file
	rewind (hotResidues);
	nLines = 0;

	while (fgets (lineString, 1000, hotResidues) != NULL)
	{
		sscanf (lineString, "%s\n", &hotResNames[nLines].resName);
		nLines++;
	}

	return hotResNames;
}

HOT_INTERACTIONS *readHotInteractionResNames (FILE *hotInteractions, HOT_INTERACTIONS *hotInteractionResNames, int *nHotInteractionResNames)
{
	// Reading the number of lines in input file
	int nLines = 0;
	char lineString[1000];
	
	while (fgets (lineString, 1000, hotInteractions) != NULL) {
		nLines++; }

	hotInteractionResNames = (HOT_INTERACTIONS *) malloc (nLines * sizeof (HOT_INTERACTIONS));
	(*nHotInteractionResNames) = nLines;

	// Storing the values from input file
	rewind (hotInteractions);
	nLines = 0;

	while (fgets (lineString, 1000, hotInteractions) != NULL) {
		sscanf (lineString, "%s %s\n", &hotInteractionResNames[nLines].resName1, &hotInteractionResNames[nLines].resName2);
		nLines++; }

	return hotInteractionResNames;
}

int main(int argc, char const *argv[])
{
	FILE *ffBondedITP, *ffNonbondedITP, *topolTopITP, *ffBondedITP_output, *ffNonbondedITP_output, *topolTopITP_output, *hotResidues, *hotInteractions;
	ffBondedITP = fopen (argv[1], "r");
	ffNonbondedITP = fopen (argv[2], "r");
	topolTopITP = fopen (argv[3], "r");
	ffBondedITP_output = fopen (argv[4], "w");
	ffNonbondedITP_output = fopen (argv[5], "w");
	topolTopITP_output = fopen (argv[6], "w");
	float lambda = atof (argv[7]);
	hotResidues = fopen (argv[8], "r");
	hotInteractions = fopen (argv[9], "r");

	// Structs to store information from topology file
	TOPOLOGY_BOOL topCurrentPosition;
	TOPOLOGY_ATOMS *inputAtoms;
	TOPOLOGY_BONDS *inputBonds;
	TOPOLOGY_PAIRS *inputPairs;
	TOPOLOGY_ANGLES *inputAngles;
	TOPOLOGY_DIHEDRALS *inputDihedrals;
	POSITIONAL_RESTRAINTS *inputRestraints;
	TOPOLOGY_SYSTEM *inputSystem;
	TOPOLOGY_MOLECULE *inputMolecule;

	// Structs to store bonded interactions
	BONDED_DEFINES *inputBondDefines, *inputAngleDefines, *inputProperDihedralDefines, *inputImproperDihedralDefines;
	int nAtoms, nBondDefines, nAngleDefines, nProperDihedralDefines, nImproperDihedralDefines;

	// Structs to store non-bonded interactions
	NONBONDED_ATOMTYPES *inputNonbondedAtomtypes;
	NONBONDED_PARAMS *inputNonbondedParams;
	int nNonbondedAtomtypes, nNonbondedParams;

	inputAtoms = readTopAtoms (topolTopITP, topCurrentPosition, inputAtoms, &nAtoms, topolTopITP_output, lambda);

	readBondedITP (ffBondedITP, topCurrentPosition, &inputBondDefines, &inputAngleDefines, &inputProperDihedralDefines, &inputImproperDihedralDefines, &nBondDefines, &nAngleDefines, &nProperDihedralDefines, &nImproperDihedralDefines, ffBondedITP_output, lambda);

	// Reading hot residues and hot interactions for Hamiltonian modification
	HOT_RESIDUES *hotResNames; int nHotResNames;
	HOT_INTERACTIONS *hotInteractionResNames; int nHotInteractionResNames;

	hotResNames = readHotResNames (hotResidues, hotResNames, &nHotResNames);
	hotInteractionResNames = readHotInteractionResNames (hotInteractions, hotInteractionResNames, &nHotInteractionResNames);

	readNonbondedITP (ffNonbondedITP, topCurrentPosition, &inputNonbondedAtomtypes, &inputNonbondedParams, &nNonbondedAtomtypes, &nNonbondedParams, ffNonbondedITP_output, lambda, hotResNames, hotInteractionResNames);

	free (inputBondDefines);
	free (inputAngleDefines);
	free (inputProperDihedralDefines);
	free (inputImproperDihedralDefines);

	free (inputAtoms);
	// free (inputBonds);
	// free (inputPairs);
	// free (inputAngles);
	// free (inputDihedrals);
	// free (inputRestraints);
	// free (inputSystem);
	// free (inputMolecule);

	fclose (ffBondedITP);
	fclose (ffNonbondedITP);
	fclose (topolTopITP);
	fclose (ffBondedITP_output);
	fclose (ffNonbondedITP_output);
	fclose (topolTopITP_output);
	return 0;
}
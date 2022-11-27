#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

typedef struct residueFraction
{
	int number;
	float fraction;
} RESIDUE_FRACTION;

float computeFraction (char lineString[], const char *identifier)
{
	float fraction = 0;
	float nTotal = 0, nFound = 0;

	// printf("length: %d\n",strlen (lineString));

	for (int i = 0; i < strlen (lineString); ++i)
	{
		if (lineString[i] != '"' && lineString[i] != ',') 
		{
			nTotal++;

			if (lineString[i] == identifier[0]) 
			{
				nFound++;
			}
		}
	}

	fraction = (float) (nFound) / (float) (nTotal - 1);
	// printf("%f/%f = %f\n", nFound, nTotal, fraction);

	return fraction;
}

int main(int argc, char const *argv[])
{
	if (argc != 5)
	{
		printf("ERROR: Incorrect args passed! The program will terminate now.\n\nREQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n{~} argv[0] = ./plotSecondaryStructures\n{~} argv[1] = Input *.xpm file name\n{~} argv[2] = Residue identifier\n{~} argv[3] = output file name\n{~} argv[4] = Number of residues in input *.xpm file.\n\nNOTE: Use escape characters at appropriate places.\n\n");
		exit (1);
	}

	FILE *input, *output;
	input = fopen (argv[1], "r");
	output = fopen (argv[3], "w");

	char lineString[10000];
	char identifierChar[5];
	snprintf (identifierChar, 1, "%s", argv[2]);

	int isResidueLine = 0, currentResidueN = 0;

	int totalResidues = atoi (argv[4]);

	RESIDUE_FRACTION *res;
	res = (RESIDUE_FRACTION *) malloc (totalResidues * sizeof (RESIDUE_FRACTION));

	while (fgets (lineString, 10000, input) != NULL)
	{
		if (strstr (lineString, "\"") && (isResidueLine == 1) && strstr (lineString, "\",") || 
			strstr (lineString, "\"") && (isResidueLine == 1) && strstr (lineString, "\"") && currentResidueN >= 35)
		{
			res[currentResidueN].fraction = computeFraction (lineString, argv[2]);
			currentResidueN++;
		}

		if (strstr (lineString, "\"") && strstr (lineString, "\",")) {
			isResidueLine = 1; }
	}

	int N = 1;
	for (int i = (currentResidueN - 1); i >= 0; --i)
	{
		fprintf(output, "%d %f\n", N, res[i].fraction);
		N++;
	}

	fclose (input);
	fclose (output);
	return 0;
}

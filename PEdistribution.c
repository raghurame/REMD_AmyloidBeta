#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

typedef struct peData
{
	float time, potential, kinetic, totalEnergy, temperature;
} PE_DATA;

typedef struct PEdistribution
{
	float potlo, pothi;
	int count;
} PE_DISTRIBUTION;

int countDatalines (FILE *input)
{
	rewind (input);
	int nDatalines = 0;
	char lineString[3000];

	while (fgets (lineString, 3000, input) != NULL)
	{
		if (lineString[0] != '@' && lineString[0] != '#')
		{
			nDatalines++;
		}
	}

	rewind (input);
	return nDatalines;
}

PE_DATA *readInputData (PE_DATA *inputData, FILE *input)
{
	char lineString[3000];
	int arrayCounter = 0;

	while (fgets (lineString, 3000, input) != NULL)
	{
		if (lineString[0] != '@' && lineString[0] != '#')
		{
			sscanf (lineString, "%f %f %f %f %f\n", 
				&inputData[arrayCounter].time, 
				&inputData[arrayCounter].potential, 
				&inputData[arrayCounter].kinetic, 
				&inputData[arrayCounter].totalEnergy, 
				&inputData[arrayCounter].temperature);
			arrayCounter++;
		}
	}

	return inputData;
}

void findMinMaxPotential (PE_DATA *inputData, int nDatalines, float *minPotential, float *maxPotential, float minTimerange, float maxTimerange)
{
	int firstLine = 1;

	for (int i = 0; i < nDatalines; ++i)
	{
		if ((inputData[i].time >= minTimerange) && (inputData[i].time <= maxTimerange))
		{
			if (firstLine == 1) {
				(*minPotential) = inputData[i].potential;
				(*maxPotential) = inputData[i].potential;
				firstLine = 0; }
			else
			{
				if (inputData[i].potential < (*minPotential)) {
					(*minPotential) = inputData[i].potential; }
				if (inputData[i].potential > (*maxPotential)) {
					(*maxPotential) = inputData[i].potential; }
			}			
		}
	}
}

void selectTimerange (PE_DATA *inputData, int nDatalines, float *minTimerange, float *maxTimerange)
{
	float beginningTime, endTime;
	for (int i = 0; i < nDatalines; ++i)
	{
		if (i == 0) {
			beginningTime = inputData[i].time; }

		endTime = inputData[i].time;
	}

	printf("The time range present in the input file: %f to %f\n", beginningTime, endTime);
	printf("Enter min time in the input range: "); scanf ("%f", &(*minTimerange)); printf("\n");
	printf("Enter max time in the input range: "); scanf ("%f", &(*maxTimerange)); printf("\n");
}

PE_DISTRIBUTION *computeDistribution (PE_DISTRIBUTION *potDistribution, PE_DATA *inputData, int nDatalines, float minPotential, float maxPotential, float potentialBinWidth, int nPotentialBins)
{
	float currentBinlo, currentBinhi;

	for (int i = 0; i < nPotentialBins; ++i)
	{
		currentBinlo = minPotential + (potentialBinWidth * i);
		currentBinhi = currentBinlo + potentialBinWidth;

		potDistribution[i].potlo = currentBinlo;
		potDistribution[i].pothi = currentBinhi;

		for (int j = 0; j < nDatalines; ++j)
		{
			if (inputData[j].potential >= currentBinlo && inputData[j].potential < currentBinhi)
			{
				potDistribution[i].count++;
			}
		}
	}

	return potDistribution;
}

void printPotDistribution (PE_DISTRIBUTION *potDistribution, int nPotentialBins, FILE *output)
{
	for (int i = 0; i < nPotentialBins; ++i)
	{
		fprintf(output, "%d %f %f %d\n", i + 1, potDistribution[i].potlo, potDistribution[i].pothi, potDistribution[i].count);
	}
}

void printAvgSDTemperature (PE_DATA *inputData, int nDatalines, float minTimerange, float maxTimerange, const char *inputFilename)
{
	float avgTemp = 0, sdTemp = 0, nDenom = 0;

	for (int i = 0; i < nDatalines; ++i)
	{
		if ((inputData[i].time >= minTimerange) && (inputData[i].time <= maxTimerange))
		{
			avgTemp += inputData[i].temperature;
			nDenom++;
		}
	}

	avgTemp /= nDenom;

	for (int i = 0; i < nDatalines; ++i)
	{
		if ((inputData[i].time >= minTimerange) && (inputData[i].time <= maxTimerange))
		{
			sdTemp += pow ((inputData[i].temperature - avgTemp), 2);
		}
	}

	sdTemp = sqrt (sdTemp / nDenom);

	printf("Avg. temperature: %f\nStdev. temperature: %f\n\n", avgTemp, sdTemp);

	FILE *outputTempStats;
	char *tempStatsFilename;
	tempStatsFilename = (char *) malloc (50 * sizeof (char));
	snprintf (tempStatsFilename, 50, "%s.temp", inputFilename);
	outputTempStats = fopen (tempStatsFilename, "w");

	fprintf(outputTempStats, "avg/stdev temperature: %f %f\n\n", avgTemp, sdTemp);

	fclose (outputTempStats);
	free (tempStatsFilename);
}

int main(int argc, char const *argv[])
{
	if (argc != 2)
	{
		printf("\nREQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n[~] argv[0] = ./program\n[~] argv[1] = input xvg file name\n\nEXPECTED COLUMNS:\n~~~~~~~~~~~~~~~~~\n1. time, 2. potential energy, 3. kinetic energy, 4. total energy, 5. temperature\n\n");
		exit (1);
	}

	FILE *input, *output;
	input = fopen (argv[1], "r");

	char *outputString;
	outputString = (char *) malloc (50 * sizeof (char));
	snprintf (outputString, 50, "%s.output", argv[1]);

	output = fopen (outputString, "w");

	int nDatalines = countDatalines (input);

	PE_DATA *inputData;
	printf("\nAllocating %d mem (PE_DATA) to store PE input data.\n\n", nDatalines);
	inputData = (PE_DATA *) malloc (nDatalines * sizeof (PE_DATA));

	inputData = readInputData (inputData, input);

	float minTimerange, maxTimerange;
	selectTimerange (inputData, nDatalines, &minTimerange, &maxTimerange);

	float minPotential, maxPotential;
	findMinMaxPotential (inputData, nDatalines, &minPotential, &maxPotential, minTimerange, maxTimerange);

	printf("Min potential in selected range: %f\nMax potential in selected range: %f\n\n", minPotential, maxPotential);
	printf("The range will be divided into 20 bins...\n");
	float potentialRange = maxPotential - minPotential;
	float potentialBinWidth = potentialRange / (float) 20;
	printf("Bin length (pot.): %f\n\n", potentialBinWidth);

	PE_DISTRIBUTION *potDistribution;
	// Allocating extra memory
	potDistribution = (PE_DISTRIBUTION *) malloc ((20 + 2) * sizeof (PE_DISTRIBUTION));

	potDistribution = computeDistribution (potDistribution, inputData, nDatalines, minPotential, maxPotential, potentialBinWidth, (int) 20 + 2);
	printPotDistribution (potDistribution, (int) 20 + 2, output);

	// Calculating avg. and SD of temperature within the selected range
	printAvgSDTemperature (inputData, nDatalines, minTimerange, maxTimerange, argv[1]);

	fclose (input);
	fclose (output);
	return 0;
}
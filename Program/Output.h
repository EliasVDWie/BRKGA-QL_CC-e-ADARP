/************************************************************************************
									IO Functions
*************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Data.h"
#include "Read.h"

void WriteSolutionScreen(TSol s, int n, float timeBest, float timeTotal, char instance[], int nbVehicles)
{
	printf("\n\n\n Instance: %s \nsol: ", instance);
	for (int i = 0; i < nbVehicles; i++) {
		printf("Vehicle %d\n", i + 1);
		printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "stop", "T", "w", "ET", "LT", "arr", "dep", "MaxRT", "C", "B");
		printf("%10d %10lf %10lf %10lf %10lf %10lf %10lf %10s %10d %10lf\n", s.sol[i][0].stop, s.sol[i][0].T, s.sol[i][0].w, s.sol[i][0].ET, s.sol[i][0].LT, node[s.sol[i][0].stop].arr, node[s.sol[i][0].stop].dep, "NA", s.sol[i][0].C, s.sol[i][0].B);
		for (int j = 1; j < s.sol[i].size(); j++) {
			if (s.sol[i][j].stop < nbUsers) {
				printf("%10d %10lf %10lf %10lf %10lf %10lf %10lf %10d %10d %10lf\n", s.sol[i][j].stop, s.sol[i][j].T , s.sol[i][j].w, s.sol[i][j].ET, s.sol[i][j].LT, node[s.sol[i][j].stop].arr, node[s.sol[i][j].stop].dep, maxRideTimes[s.sol[i][j].stop], s.sol[i][j].C, s.sol[i][j].B);
			}
			else {
				printf("%10d %10lf %10lf %10lf %10lf %10lf %10lf %10s %10d %10lf\n", s.sol[i][j].stop, s.sol[i][j].T, s.sol[i][j].w, s.sol[i][j].ET, s.sol[i][j].LT, node[s.sol[i][j].stop].arr, node[s.sol[i][j].stop].dep, "NA", s.sol[i][j].C, s.sol[i][j].B);
			}
		}
	}
	printf("\nLocal search applied: %d", s.flag);
	printf("\nDecoder: %.2lf", s.vec[n].rk);
	printf("\nScheduler: %.2lf", s.vec[n - 1].rk);
	printf("\nCharging priority: %.2lf", s.vec[n - 2].rk);
	printf("\nfo: %.5lf",s.fo);
	printf("\nTotal time: %.3f",timeTotal);
	printf("\nBest time: %.3f\n\n",timeBest);
	
    //printf("\nH1: %d (%d) H2: %d (%d) H3: %d (%d) H4: %d (%d)\n\n", H1, SH1, H2, SH2, H3, SH3, H4, SH4);
}

void WriteSolution(TSol s, int n, float timeBest, float timeTotal, char instance[])
{
	FILE *arquivo;
    	arquivo = fopen("Results/Solutions.txt","a");

	if (!arquivo)
	{
		printf("\n\nFile not found Solutions.txt!!!");
		getchar();
		exit(1);
	}

	fprintf(arquivo, "\n\n\n Instance: %s \nsol: ", instance);
	for (int i = 0; i < nbVehicles; i++) {
		fprintf(arquivo, "Vehicle %d\n", i + 1);
		fprintf(arquivo, "%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "stop", "T", "w", "ET", "LT", "arr", "dep", "MaxRT", "C", "B");
		fprintf(arquivo, "%10d %10lf %10lf %10lf %10lf %10lf %10lf %10s %10d %10lf\n", s.sol[i][0].stop, s.sol[i][0].T, s.sol[i][0].w, s.sol[i][0].ET, s.sol[i][0].LT, node[s.sol[i][0].stop].arr, node[s.sol[i][0].stop].dep, "NA", s.sol[i][0].C, s.sol[i][0].B);
		for (int j = 1; j < s.sol[i].size(); j++) {
			if (s.sol[i][j].stop < nbUsers) {
				fprintf(arquivo, "%10d %10lf %10lf %10lf %10lf %10lf %10lf %10d %10d %10lf\n", s.sol[i][j].stop, s.sol[i][j].T, s.sol[i][j].w, s.sol[i][j].ET, s.sol[i][j].LT, node[s.sol[i][j].stop].arr, node[s.sol[i][j].stop].dep, maxRideTimes[s.sol[i][j].stop], s.sol[i][j].C, s.sol[i][j].B);
			}
			else {
				fprintf(arquivo, "%10d %10lf %10lf %10lf %10lf %10lf %10lf %10s %10d %10lf\n", s.sol[i][j].stop, s.sol[i][j].T, s.sol[i][j].w, s.sol[i][j].ET, s.sol[i][j].LT, node[s.sol[i][j].stop].arr, node[s.sol[i][j].stop].dep, "NA", s.sol[i][j].C, s.sol[i][j].B);
			}
		}
	}
	fprintf(arquivo, "\nLocal search applied: %d", s.flag);
	fprintf(arquivo, "\nDecoder: %.2lf", s.vec[n].rk);
	fprintf(arquivo, "\nScheduler: %.2lf", s.vec[n - 1].rk);
	fprintf(arquivo, "\nCharging priority: %.2lf", s.vec[n - 2].rk);
	fprintf(arquivo, "\nfo: %.5lf", s.fo);
	fprintf(arquivo, "\nTotal time: %.3f", timeTotal);
	fprintf(arquivo, "\nBest time: %.3f\n\n", timeBest);

	fclose(arquivo);
}

void WriteResults(double fo, double foAverage, std::vector <double> fos, std::vector <bool> lss, std::vector <float> decs, std::vector <float> scheds, std::vector <float> cps, float timeBest, float timeTotal, char instance[], int ls, double TRT, double ERT, double CC)
{
	FILE *arquivo;
    	arquivo = fopen("Results/Results.csv","a");

	if (!arquivo)
	{
		printf("\n\nFile not found Results.xls!!!");
		getchar();
		exit(1);
	}

	fprintf(arquivo,"\n%s", instance);
	fprintf(arquivo, "\t%d", ls);
    	fprintf(arquivo,"\t%d", (int)fos.size());
    	for (int i=0; i<fos.size(); i++)
        	fprintf(arquivo,"\t%lf \t%s \t%f \t%f \t%f", fos[i], lss[i]?"Yes":"No", decs[i], scheds[i], cps[i]);
	fprintf(arquivo,"\t%lf", fo);
	fprintf(arquivo, "\t%lf", TRT);
	fprintf(arquivo, "\t%lf", ERT);
	fprintf(arquivo, "\t%lf", CC);
	fprintf(arquivo,"\t%lf", foAverage);
	fprintf(arquivo,"\t%.3f", timeBest);
	fprintf(arquivo,"\t%.3f", timeTotal);

	fclose(arquivo);
}

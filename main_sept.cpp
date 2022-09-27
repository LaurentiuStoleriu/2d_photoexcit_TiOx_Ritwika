#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "alglibmisc.h"
#include <algorithm>
#include <vector>
#include <random>

using namespace alglib;

const int nPart = 81; // 40000;
const int nSide = 9; // 200;
const int nMaxNeigh = 4;

const int nMaxSteps = 1000000;

const double radius = 0.2;
const double radiusLS = radius;
const double radiusHS = 0.22;
const double L = 0.6;
double depth = 0.0;
int nH, nL;

const double tempLimDown = 300.0;
const double tempLimUp = 470.0;
const double tempExcitation = 2600.0;

// const double coefExoTerm = 10;	// deg. increase temp of each neighbour
// const double coefTerm = 0.005;    //% of temperature difference exchanged at each step
const double coefTermExt = 0.001;// 0.0005;// 0.000005;
const double coefTerm_lambda = 0.01;
const double coefTerm_beta = 0.01;
const double coefTerm_s = 0.01;

double deltaQ;
double Cp_beta = 154.25;
double Cp_lambda = 161.82;
double Cp_air = 29.0;

struct sReadData
{
	double x, y, z, r, T;
	int idxCryst;
};

struct sReadData Medium[nPart];
int neighbours[nPart][nMaxNeigh];
int noOfNeighbours[nPart];
double tempAtBegin[nPart];

//char sysFile[500]        = "/home/ritwika/data/1.hc4250_Oct'20/results/systems/200x200_r02_d06.dat";
//char sysFileExcited[500] = "/home/ritwika/data/1.hc4250_Oct'20/results/200x200_r02_d06_excit.dat";
//char rezFile[500]        = "/home/ritwika/data/1.hc4250_Oct'20/results/200x200_T2600_f1_H10_c0p01.dat";
char sysFile[500] = "E:\\Stoleriu\\C\\special\\3d\\generare\\2022\\TiOX\\9x9_r02_d06.dat";
char sysFileExcited[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\9x9_r02_d06_excit.dat";
char rezFile[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\9x9_T2600_f1_H10_c0p01.dat";
char pathUCD[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\9x9_T2600_f1_H10_c0p01";

/////////////////////////////////////////// Prototypes
void initialization(void);
void alglibFunctionNeighbours(void);
void photoExcitation(void);
double fInvers(double x);
void tempExchange(void);
void saveUCD(char* fisSaveVis);
//////////////////////////////////////////////////////

int main()
{
	char fisSaveVis[500];

	initialization();

	//photoExcitation();
	for (int i = 0; i < nPart; i++)
	{
		if ( (i == 39) || (i == 41) )
		{
			Medium[i].r = radiusHS;
			Medium[i].T = tempExcitation;
			nH++; nL--;
		}
		else
		{
			Medium[i].r = radiusLS;
			Medium[i].T = tempLimDown;
		}
	}

	FILE* fp;
	fp = fopen(sysFileExcited, "w");
	for (int i = 0; i < nPart; i++)
	{
		fprintf(fp, "%20.15lf  %20.15lf  %20.15lf  %20.15lf   %20.15lf\n", Medium[i].x, Medium[i].y, Medium[i].z, Medium[i].r, Medium[i].T);
	}
	fclose(fp);

	//*** RELAXATION
	double timeInit = 0.0;
	double stepTime = 0.1;
	double sysTime = timeInit;
	int stepCount = 0;

	FILE* frez;
	frez = fopen(rezFile, "w");

	while ((stepCount < nMaxSteps) && (nH > 0))
	{
		sysTime += stepTime;

		tempExchange();

		for (int i = 0; i < nPart; i++)
		{
			tempAtBegin[i] = Medium[i].T;
		}

		for (int i = 0; i < nPart; i++)
		{
			if ((Medium[i].r > 1.05 * radiusLS) && (tempAtBegin[i] <= tempLimUp))
			{
				Medium[i].r = radiusLS;

				//				Medium[i].T += CoefExoTerm * nMaxNeigh;			// exothermic H-to-L either by heating the particle
				// 				deltaQ = coefExoTerm * 1.0;						//                   or its neighbours
				deltaQ = 10.0; 									//trying with real values
				for (int j = 0; j < noOfNeighbours[i]; j++)
				{
					Medium[neighbours[i][j]].T += (Medium[neighbours[i][j]].r > radiusLS * 1.01) ? deltaQ / Cp_lambda : deltaQ / Cp_beta;
					//	Medium[i].T -= deltaQ/Cp_lambda;
				}

				nH--; nL++;
			}
			else
			{
				if ((Medium[i].r < 1.05 * radiusLS) && (tempAtBegin[i] >= tempLimUp))
				{
					Medium[i].r = radiusHS;
					// 					deltaQ = coefExoTerm * 1.0;
					// 					Medium[neighbours[i][j]].T -= (Medium[neighbours[i][j]].r > radiusLS*1.01) ? deltaQ/Cp_lambda : deltaQ/Cp_beta;
					deltaQ = 10.0; 	//trying with real values
					for (int j = 0; j < noOfNeighbours[i]; j++)			//                   or its neighbours
					{
						Medium[neighbours[i][j]].T -= (Medium[neighbours[i][j]].r > radiusLS * 1.01) ? deltaQ / Cp_lambda : deltaQ / Cp_beta;
						//Medium[i].T += deltaQ/Cp_beta;
					}

					nH++; nL--;
				}
			}
		}

		if (!((int)(sysTime / stepTime) % 10))
		{
			printf("Time %5.2lf \t Temp %5.2lf \t HS %d \n", sysTime, Medium[0].T, nH);
			sprintf(fisSaveVis, "%s_%07d.inp", pathUCD, (int)(sysTime / stepTime));
			saveUCD(fisSaveVis);
		}
		fprintf(frez, "%20.15lf   %d\n", sysTime, nH);
	}
	fclose(frez);

	return 0;
}

//////////////////////////////////////////////////////////////////////////////

void initialization(void)
{
	FILE* fp;
	fp = fopen(sysFile, "r");
	nH = 0;
	nL = 0;

	for (int i = 0; i < nPart; i++)
	{
		fscanf(fp, "%lf  %lf  %lf  %lf  %d\n", &Medium[i].x, &Medium[i].y, &Medium[i].z, &Medium[i].r, &Medium[i].idxCryst);
		Medium[i].T = tempLimDown;
		nL++;
		if (Medium[i].x > depth)	// to remember maximum value for x
			depth = Medium[i].x;
	}
	fclose(fp);

	alglibFunctionNeighbours();
}

//////////////////////////////////////////////////////////////////////////////

void alglibFunctionNeighbours(void)
{
	int	neighboursMax = 0, neighboursMed = 0;
	int i, j, local_index;
	double distance;
	real_2d_array a;

	a.setlength(nPart, 3);
	for (i = 0; i < nPart; i++)
	{
		a(i, 0) = Medium[i].x;
		a(i, 1) = Medium[i].y;
		a(i, 2) = Medium[i].z;
	}

	integer_1d_array tags;
	tags.setlength(nPart);
	for (int i = 0; i < nPart; i++)
	{
		tags(i) = i;
	}

	ae_int_t nx = 3;
	ae_int_t ny = 0;
	ae_int_t normtype = 2;

	kdtree kdt;
	//kdtreebuild(a, nx, ny, normtype, kdt);
	kdtreebuildtagged(a, tags, nx, ny, normtype, kdt);

	real_1d_array x;
	x.setlength(3);
	real_2d_array r = "[[]]";

	integer_1d_array indexes;

	for (i = 0; i < nPart; i++)
	{
		x(0) = Medium[i].x;
		x(1) = Medium[i].y;
		x(2) = Medium[i].z;

		ae_int_t k;

		distance = 1.1 * (2.0 * radius + L);

		k = kdtreequeryrnn(kdt, x, distance, false);

		noOfNeighbours[i] = (int)k;

		// 		if (noOfNeighbours[i] + 1 > n_max_vec - 1)
		// 		{
		// 			printf("\n\n TOO MANY NEIGHBOURS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n\n");
		// 			break;
		// 		}
		// 		else
		{
			kdtreequeryresultstags(kdt, indexes);

			//SORT FOR TYPE
			std::vector<int> myvector(noOfNeighbours[i]);
			for (j = 0; j < noOfNeighbours[i]; j++)
				myvector[j] = (int)indexes(j);

			std::sort(myvector.begin(), myvector.end());

			for (j = 0; j < noOfNeighbours[i]; j++)
			{
				local_index = myvector[j];// (int)indexes(j);

				neighbours[i][j] = local_index;  //<<<---- asta e vecin
			}
			//n_vec++;
		}

		//printf("particula: %d  cu  %d  vecini \n", i, neighbours[i]);
		neighboursMed += noOfNeighbours[i];
		if (noOfNeighbours[i] > neighboursMax)	neighboursMax = noOfNeighbours[i];
	}
	printf("Maximum %d neighbours, an average of %f neighbours\n", neighboursMax, (double)neighboursMed / nPart);
	//getchar();
}

//////////////////////////////////////////////////////////////////////////////

void photoExcitation(void)
{
	double valueToCheck;

	std::random_device rd;
	//std::mt19937_64 gen(rd());    // radnom seed
	std::mt19937_64 gen(1);   // fixed seed
	std::uniform_real_distribution<double> rand_dis(0.0, 1.0);  // use with rand_dis(gen)

	// for (int i = 0; i < 10; i++)
	// {
	//     printf("%lf\n", rand_dis(gen));
	// }

	for (int fluency = 0; fluency < 1; fluency++)
	{
		for (int i = 0; i < nPart; i++)
		{
			valueToCheck = fInvers((Medium[i].x / depth) + 0.01) + 0.01;

			if (Medium[i].r > 1.001 * radiusLS)
				continue;//if already HS go to next - conseq on fluency?

			if (valueToCheck > rand_dis(gen))
			{
				Medium[i].r = radiusHS;
				Medium[i].T = tempExcitation;
				nH++; nL--;
			}
			else
			{
				Medium[i].r = radiusLS;
				Medium[i].T = tempLimDown;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////

double fInvers(double x)
{
	//return( log(1.0 + x * (exp(depth) - 1.0)));//log(1.0 + x * (exp(depth) - 1.0)) scalat la depth = 3.0;
	//return( log(1.0 + (x * 2.0 / depth) * 6.3890560989306502272) / 2.6230812603996638992);
	return(-log(x) / 10.0);
}

//////////////////////////////////////////////////////////////////////////

void tempExchange(void)
{
	int i, j, n;
	double Q = 0.0;

	for (int i = 0; i < nPart; i++)
	{
		tempAtBegin[i] = Medium[i].T;
	}

	for (i = 0; i < nPart; i++)
	{
		if (noOfNeighbours[i] < nMaxNeigh)
		{
			deltaQ = (tempAtBegin[i] - tempLimDown) * (nMaxNeigh - noOfNeighbours[i]) * coefTermExt;
			Medium[i].T -= (Medium[i].r > radiusLS * 1.01) ? deltaQ / Cp_lambda : deltaQ / Cp_beta;
			for (j = 0; j < noOfNeighbours[i]; j++)
			{
				n = neighbours[i][j];

				if ((Medium[i].r > radiusLS * 1.01) && (Medium[n].r > radiusLS * 1.01))
				{
					Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_lambda;
				}
				else
				{
					if ((Medium[i].r > radiusLS * 1.01) && (Medium[n].r < radiusHS * 1.01))
					{
						Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_s;
					}
					else
					{
						if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r > radiusLS * 1.01))
						{
							Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_s;
						}
						else
						{
							//if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r < radiusHS * 1.01))
							//{
								Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_beta;
							//}
						}
					}
				}
				Medium[i].T -= (Medium[i].r > radiusLS * 1.01) ? Q / Cp_lambda : Q / Cp_beta;
				Medium[n].T += (Medium[n].r > radiusLS * 1.01) ? Q / Cp_lambda : Q / Cp_beta;
			}
		}
		else
		{
			for (j = 0; j < noOfNeighbours[i]; j++)
			{
				n = neighbours[i][j];
				if ((Medium[i].r > radiusLS * 1.01) && (Medium[n].r > radiusLS * 1.01))
				{
					Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_lambda;
				}
				else
				{
					if ((Medium[i].r > radiusLS * 1.01) && (Medium[n].r < radiusHS * 1.01))
					{
						Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_s;
					}
					else
					{
						if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r > radiusLS * 1.01))
						{
							Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_s;
						}
						else
						{
							//if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r < radiusHS * 1.01))
							//{
								Q = (tempAtBegin[i] - tempAtBegin[n]) * coefTerm_beta;
							//}
						}
					}
				}
				//Q = ((Medium[i].r > radiusLS*1.01) && (Medium[n].r < radiusHS*1.01)) ? (Medium[i].T - Medium[n].T)*coefTerm_s : (Medium[i].T - Medium[n].T)*coefTerm_lambda;

				Medium[i].T -= (Medium[i].r > radiusLS * 1.01) ? Q / Cp_lambda : Q / Cp_beta;
				Medium[n].T += (Medium[n].r > radiusLS * 1.01) ? Q / Cp_lambda : Q / Cp_beta;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////

void saveUCD(char* fisSaveVis)
{
	int p, i, j, v1, v2;
	int count = 0;
	int* count_switched;
	double d1;

	count_switched = (int*)calloc(nPart, sizeof(int));

	for (p = 0; p < nPart; p++)
	{
		count_switched[p] = 0;

		for (i = 0; i < noOfNeighbours[p]; i++)
		{
			if (Medium[neighbours[p][i]].r > radiusLS * 1.01)
			{
				count_switched[p]++;
			}
		}

		for (i = 0; i < noOfNeighbours[p] - 1; i++)
		{
			for (j = i + 1; j < noOfNeighbours[p]; j++)
			{
				v1 = neighbours[p][i];
				v2 = neighbours[p][j];

				d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

				if ((d1 < 1.9)) //to avoid two neighbours both on Ox or Oy
				{
					count++;
				}
			}
		}
	}

	//char fis_save_vis[500];
	//sprintf(fis_save_vis, "%s_ucd_%06d.inp", file, (int)timp);

	FILE* fpout;
	fpout = fopen(fisSaveVis, "w");

	fprintf(fpout, "%d %d 1 0 0\n", nPart, count);
	printf("SAVING UCD %d %d\n", nPart, count);

	for (i = 0; i < nPart; i++)
	{
		fprintf(fpout, "%d %f %f %f\n", i + 1, Medium[i].x, Medium[i].y, Medium[i].z);
	}

	count = 0;
	for (p = 0; p < nPart; p++)
	{
		for (i = 0; i < noOfNeighbours[p] - 1; i++)
		{
			for (j = i + 1; j < noOfNeighbours[p]; j++)
			{
				v1 = neighbours[p][i];
				v2 = neighbours[p][j];

				d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

				if ((d1 < 1.9)) //to avoid two neighbours both on Ox or Oy
				{
					count++;
					fprintf(fpout, "%d 1 tri  %d  %d  %d \n", count, p + 1, v1 + 1, v2 + 1);
				}
			}
		}
	}

	fprintf(fpout, "3 1 1 1\n");
	fprintf(fpout, "radius, nm\n");
	fprintf(fpout, "crystallite, au\n");
	fprintf(fpout, "Temperature, K\n");

	for (i = 0; i < nPart; i++)
	{
		fprintf(fpout, "%d %lf %d %f\n", i + 1, Medium[i].r, Medium[i].idxCryst, Medium[i].T);
	}

	fclose(fpout);

	free(count_switched);
}

//////////////////////////////////////////////////////////////////////////
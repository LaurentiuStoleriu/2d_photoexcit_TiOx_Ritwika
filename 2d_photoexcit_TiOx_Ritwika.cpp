#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "alglibmisc.h"
#include <algorithm>
#include <vector>
#include <random>

using namespace alglib;

const int nPart = 39800;
const int nMaxNeigh = 6;

const char   sysFile[500] = "E:\\Stoleriu\\C\\special\\3d\\generare\\2022\\Elastic\\100x100_RektHex_L06_LS.dat";
const char   volFile[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\100x100_RektHex_Sol_TRY.dat";

const double radius = 0.20;
const double L = 0.60;
const double radiusLS = 0.20;
const double radiusHS = 0.22;
	  double depth = 0.0;

const double tempLimDown = 200.0;
const double tempLimUp = 350.0;
const double tempExcitation = 900.0;

const double coefExoTerm = 1.0;	// deg. crestere temp a fiecari vecin
const double coefTerm = 0.5;    //% din diferenta de temperaturi ce se schimba per pas
const double coefTermExt = 0.01;

const int	 nMaxSteps = 1000000;


typedef struct
{
	public: double x, y, z, radius, theta, k;
}sReadData;

sReadData	Medium[nPart];
int			noOfNeighbours[nPart];
int			neighbours[nPart][nMaxNeigh];
double		temp[nPart];
double		tempAtBegin[nPart];

//////////////////////////////////////////////////////////////////////////
int initialization(void);
void alglibFunctionNeighbours(void);
double fInvers(double x);
int tempExchange(void);
//////////////////////////////////////////////////////////////////////////


int main()
{
	initialization();
	
	FILE *fVol;
	fVol = fopen(volFile, "w");

	int nH = 0;
	int nL = nPart;
	int nL_Old = nPart;
	double valueToCheck;

	int i, j;

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> rand_dis(0.0, 1.0);  // use with rand_dis(gen)

	//*** PHOTOEXCITATION
	for (int fluency = 0; fluency < 1; fluency++)
	{
		for (i = 0; i < nPart; i++)
		{
			valueToCheck = fInvers((Medium[i].x / depth) + 0.01) + 0.01;

			if (Medium[i].radius > 1.001 * radiusLS)
				continue;//if already HS go to next - conseq on fluency?

			if (valueToCheck > rand_dis(gen))
			{
				temp[i] = tempExcitation;
				Medium[i].radius = radiusHS;
				nH++; nL--;
			}
			else
			{
				temp[i] = tempLimDown;
				Medium[i].radius = radiusLS;
			}
		}
	}

	double tInit = 0.0;
	double stepT = 0.1;
	double sysTime = tInit;
	int stepCount = 0;
	double aria = 0.0;


	//*** RELAXATION
// 	int wait[nPart];
//  	for (i = 0; i < nPart; i++)
// 		wait[i] = 0;
	while ( (stepCount < nMaxSteps) && (nH > 0) )
	{
		sysTime += stepT;

		tempExchange();

		for (i = 0; i < nPart; i++)
			tempAtBegin[i] = temp[i];

		for (i = 0; i < nPart; i++)
		{
			if ((Medium[i].radius > 1.05 * radiusLS) && (tempAtBegin[i] <= tempLimUp))
			{
//				wait[i] = 5;
				Medium[i].radius = radiusLS;

//				temp[i] += CoefExoTerm * 6.0;					// exothermic H-to-L either by heating the particle
				for (j = 0; j < noOfNeighbours[i]; j++)			//                   or its neighbours
				{
					temp[neighbours[i][j]] += coefExoTerm * 1.0;
				}

				nH--; nL++;
			}
			else
			{
				if ((Medium[i].radius < 1.05 * radiusLS) && (tempAtBegin[i] >= tempLimUp) /*&& (!wait[i])*/)
				{
					Medium[i].radius = radiusHS;
					nH++; nL--;
				}
			}
// 			if (wait[i] > 0)
// 				wait[i]--;
		}

		if (((stepCount <= 1000) && !(stepCount % 100)) || ((stepCount > 1000) && !(stepCount % 1000)))
		{
#ifdef grafic
			{
				count = 0;
				count_switched = (int*)calloc(n_part, sizeof(int));

				for (int p = 0; p < n_part; p++)
				{
					count_switched[p] = 0;

					for (int i = 0; i < neighbours[p]; i++)
					{
						if (Medium[Position_Coef[p][i].vecin].raza > 1.05)
						{
							count_switched[p]++;
						}
					}

					for (int i = 0; i < neighbours[p] - 1; i++)
					{
						for (int j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < L + 3.0 * radius)) //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
							}
						}
					}
				}


				char fis_save_vis[500];
				sprintf(fis_save_vis, "%s_ucd_%06d.inp", file, (int)timp);

				FILE* fpout;
				fpout = fopen(fis_save_vis, "w");

				fprintf(fpout, "%d %d 1 0 0\n", n_part, count);
				printf("SAVING UCD %d %d\n", n_part, count);

				for (int i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %f %f %f\n", i + 1, Medium[i].x, Medium[i].y, Medium[i].z);
				}

				count = 0;
				for (int p = 0; p < n_part; p++)
				{
					for (int i = 0; i < neighbours[p] - 1; i++)
					{
						for (int j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < L + 3.0 * radius))  //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
								fprintf(fpout, "%d 1 tri  %d  %d  %d \n", count, p + 1, v1 + 1, v2 + 1);
							}

						}
					}
				}

				fprintf(fpout, "6 1 1 1 1 1 1\n");
				fprintf(fpout, "raza, nm\n");
				fprintf(fpout, "phase, au\n");
				fprintf(fpout, "PLH, au\n");
				fprintf(fpout, "PHL, au\n");
				fprintf(fpout, "Temperature, K\n");
				fprintf(fpout, "pressure, au\n");

				for (int i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %lf %lf %lf %lf %lf %lf\n", i + 1, Medium[i].raza, Medium[i].k,
						((probabilitateLH[i] < 1.0) ? (probabilitateLH[i]) : (1.0)),
						((probabilitateHL[i] < 1.0) ? (probabilitateHL[i]) : (1.0)),
						T[i],
						pres[i]);
				}

				fclose(fpout);

				free(count_switched);
			}

#endif
#ifdef ThermalHisto
			doTheHisto((int)timp);
#endif
			printf("Time %5.2lf \t Temp %5.2lf \t HS %d \n", sysTime, temp[0], nH);
		}


	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////

int initialization(void)
{
	FILE* fp;
	long i;

	/// READ Medium
	fp = fopen(sysFile, "r");
	for (i = 0; i < nPart; i++)
	{
		fscanf(fp, "%lG %lG %lG %lG %lG %lG \n", &Medium[i].x, &Medium[i].z, &Medium[i].y, &Medium[i].radius, &Medium[i].theta, &Medium[i].k);
		if (Medium[i].x > depth)	// to remember maximum value for x
			depth = Medium[i].x;
	}
	fclose(fp);

	///// COMPUTE NEIGHBOURS
	alglibFunctionNeighbours();

	///// PRINT NEIGHBOURS
// 	for (i = 0; i < n_part; i++)
// 	{
// 		printf("%d: ", i);
// 		for (j = 0; j < neighbours[i]; j++)
// 		{
// 			printf(" %d ", Position_Coef[i][j].vecin);
// 		}
// 		printf("\n");
// 	}

	printf("Read %d particles \n", nPart);

	return(0);
}

//////////////////////////////////////////////////////////////////////////

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
		//k = kdtreequeryknn(kdt, x, 2, false);
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

//////////////////////////////////////////////////////////////////////////

double fInvers(double x)
{
	//return( log(1.0 + x * (exp(depth) - 1.0)));//log(1.0 + x * (exp(depth) - 1.0)) scalat la depth = 3.0;
	//return( log(1.0 + (x * 2.0 / depth) * 6.3890560989306502272) / 2.6230812603996638992); 
	return(-log(x) / 10.0);
}

//////////////////////////////////////////////////////////////////////////

int tempExchange(void)
{
	int i, j, v;
	double Q;
	for (i = 0; i < nPart; i++)
	{

		if (noOfNeighbours[i] < 5)
		{
			temp[i] -= (temp[i] - tempLimDown) * (5 - noOfNeighbours[i]) * coefTermExt;
			//T[i] = T_LIM_DWN;
		}
		else
		{
			for (j = 0; j < noOfNeighbours[i]; j++)
			{
				v = neighbours[i][j];
				Q = (temp[i] - temp[v]) * coefTerm;
				temp[i] -= Q;
				temp[v] += Q;
			}
		}

	}
	return(0);
}

//////////////////////////////////////////////////////////////////////////


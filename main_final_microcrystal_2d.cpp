#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "alglibmisc.h"
#include <algorithm>
#include <vector>
#include <random>

using namespace alglib;

    const int nPart = 40000;
    const int nSide = 200;
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
	const double coefTermExt = 0.1;// 0.0005;// 0.000005; 
	const double coefTerm_lambda = 1;
	const double coefTerm_beta = 2;
	const double coefTerm_s = 0.7;
		
	double deltaQ;
	double Cp_beta = 154.25; 
	double Cp_lambda = 161.82;
	double Cp_air = 29;



    struct sReadData
    {
        double x, y, z, r, T;
		int idxCryst;
    };

    struct sReadData Medium[nPart];
    int neighbours[nPart][nMaxNeigh];
    int noOfNeighbours[nPart];
    double tempAtBegin[nPart];

//     char sysFile[500]        = "/home/ritwika/data/1.hc4250_Oct'20/results/systems/200x200_r02_d06.dat";
//     char sysFileExcited[500] = "/home/ritwika/data/1.hc4250_Oct'20/results/200x200_r02_d06_excit.dat";
//     char rezFile[500]        = "/home/ritwika/data/1.hc4250_Oct'20/results/200x200_T2600_f0.dat";
    char sysFile[500]        = "E:\\Stoleriu\\C\\special\\3d\\generare\\2022\\TiOX\\200x200_r02_d06.dat";
    char sysFileExcited[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\200x200_r02_d06_excit.dat";
    char rezFile[500]        = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\200x200_T2600_f1_cl1_cb10.dat";


/////////////////////////////////////////// Prototypes
void initialization(void);
void alglibFunctionNeighbours(void);
void photoExcitation(void);
double fInvers(double x);
void tempExchange(void);
//////////////////////////////////////////////////////

int main()
{
    initialization();

    photoExcitation();

    FILE *fp;
    fp = fopen(sysFileExcited, "w");
    for (int i = 0; i < nPart; i++)
    {
        fprintf(fp, "%20.15lf  %20.15lf  %20.15lf  %20.15lf   %20.15lf  %d\n", Medium[i].x, Medium[i].y, Medium[i].z, Medium[i].r, Medium[i].T, Medium[i].idxCryst);
    }
    fclose(fp);

	//*** RELAXATION
    double timeInit = 0.0;
	double stepTime = 0.1;
	double sysTime = timeInit;
	int stepCount = 0;

    FILE *frez;
	frez = fopen(rezFile, "w");

    while ( (stepCount < nMaxSteps) && (nH > 0) )
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
				for (int j = 0; j < noOfNeighbours[i]; j++)			//                   or its neighbours
				{
					//deltaQ = coefExoTerm * 1.0;
					deltaQ = 12000; 	//trying with real values
					Medium[neighbours[i][j]].T += (Medium[neighbours[i][j]].r > radiusLS*1.01) ? deltaQ/Cp_lambda : deltaQ/Cp_beta;

				//	Medium[i].T -= deltaQ/Cp_lambda;
				}

				nH--; nL++;
			}
			else
			{
				if ((Medium[i].r < 1.05 * radiusLS) && (tempAtBegin[i] >= tempLimUp))
				{
					Medium[i].r = radiusHS;
					for (int j = 0; j < noOfNeighbours[i]; j++)			//                   or its neighbours
					{
						//deltaQ = coefExoTerm * 1.0;
						//Medium[neighbours[i][j]].T -= (Medium[neighbours[i][j]].r > radiusLS*1.01) ? deltaQ/Cp_lambda : deltaQ/Cp_beta;
						deltaQ = 12000 ; 	//trying with real values
						Medium[neighbours[i][j]].T -= (Medium[neighbours[i][j]].r > radiusLS*1.01) ? deltaQ/Cp_lambda : deltaQ/Cp_beta;
						//Medium[i].T += deltaQ/Cp_beta;
					}	
					
					nH++; nL--;
				}
			}
		}

		stepCount++;
		if (!(stepCount % 100))
		{
			printf("Time %5.2lf \t Temp %5.2lf \t HS %d \n", sysTime, Medium[0].T, nH);
		}
        fprintf(frez, "%20.15lf \t %d \t %5.2lf\n", sysTime, nH, Medium[0].T);
    }
    fclose(frez);


    return 0;
}

//////////////////////////////////////////////////////////////////////////////

void initialization(void)
{
    FILE *fp;
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
	//std::mt19937_64 gen(rd());    // random seed
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
	double Q;

	for (i = 0; i < nPart; i++)
	{
		if (noOfNeighbours[i] < nMaxNeigh)
		{
			deltaQ = (Medium[i].T - tempLimDown) * (nMaxNeigh - noOfNeighbours[i]) * coefTermExt;
			Medium[i].T -= (Medium[i].r > radiusLS * 1.01) ? deltaQ / Cp_lambda : deltaQ / Cp_beta;
			for (j = 0; j < noOfNeighbours[i]; j++)
			{
				n = neighbours[i][j];

				if ((Medium[i].r > radiusLS * 1.01) && (Medium[n].r > radiusLS * 1.01))
				{
					Q = (Medium[i].T - Medium[n].T) * coefTerm_lambda;
				}
				else
				{
					if ((Medium[i].r > radiusLS * 1.01) && (Medium[n].r < radiusHS * 1.01))
					{
						Q = (Medium[i].T - Medium[n].T) * coefTerm_s;
					}
					else
					{
						if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r > radiusLS * 1.01))
						{
							Q = (Medium[i].T - Medium[n].T) * coefTerm_s;
						}
						else
						{
							if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r < radiusHS * 1.01))
							{
								Q = (Medium[i].T - Medium[n].T) * coefTerm_beta;
							}
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
					Q = (Medium[i].T - Medium[n].T) * coefTerm_lambda;
				}
				else
				{
					if ((Medium[i].r > radiusLS * 1.01) && (Medium[n].r < radiusHS * 1.01))
					{
						Q = (Medium[i].T - Medium[n].T) * coefTerm_s;
					}
					else
					{
						if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r > radiusLS * 1.01))
						{
							Q = (Medium[i].T - Medium[n].T) * coefTerm_s;

						}
						else
						{
							if ((Medium[i].r < radiusHS * 1.01) && (Medium[n].r < radiusHS * 1.01))
							{
								Q = (Medium[i].T - Medium[n].T) * coefTerm_beta;
							}
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
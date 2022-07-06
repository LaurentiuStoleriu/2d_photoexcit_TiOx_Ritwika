#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "alglibmisc.h"
#include <algorithm>
#include <vector>
using namespace alglib;


const int nPart = 39800;
const int nMaxNeigh = 6;

const double radius = 0.2;
const double L = 0.6;
	  double depth = 0.0;

constexpr char sys_file[500] = "E:\\Stoleriu\\C\\special\\3d\\generare\\2022\\Elastic\\100x100_RektHex_L06_LS.dat"; // L = 0.6, LS:r = 0.2, HS:r=0.22

typedef struct
{
	public: double x, y, z, raza, theta, k;
}sReadData;

sReadData	Medium[nPart];
int			noOfNeighbours[nPart];
int			neighbours[nPart][nMaxNeigh];

//////////////////////////////////////////////////////////////////////////
int initialization(void);
void alglibFunctionNeighbours(void);
//////////////////////////////////////////////////////////////////////////


int main()
{
    initialization();
	///////////////////////////////////////////////////////////////////////


}


//////////////////////////////////////////////////////////////////////////

int initialization(void)
{
	FILE* fp;
	long i;

	/// READ Medium
	fp = fopen(sys_file, "r");
	for (i = 0; i < nPart; i++)
	{
		fscanf(fp, "%lG %lG %lG %lG %lG %lG \n", &Medium[i].x, &Medium[i].z, &Medium[i].y, &Medium[i].raza, &Medium[i].theta, &Medium[i].k);
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
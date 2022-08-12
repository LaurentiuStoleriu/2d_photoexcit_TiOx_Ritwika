#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>

const int nPart = 10000;
const int nSide = 100;
const int nSideCryst = 10;

double radius = 0.2;
double distance = 0.6;

struct sReadData
{
    double x, y, z, r;
    int idxCryst;
};

struct sReadData Medium[nPart];

int main()
{
    double stepX = distance + 2.0 * radius;
    double stepY = distance + 2.0 * radius;
    double x0 = radius;
    double y0 = radius;
    for (int i = 0; i < nSide; i++)
    {
        for (int j = 0; j < nSide; j++)
        {
            Medium[i*nSide + j].x = x0 + j * stepX;
            Medium[i*nSide + j].y = y0 + i * stepY;
            Medium[i*nSide + j].z = 0.0;
            Medium[i*nSide + j].r = radius;
            Medium[i*nSide + j].idxCryst = 1 + (i/nSideCryst)*(nSide/nSideCryst) + j/nSideCryst;
        }
        
    }

    FILE *fp;
    //fp = fopen("/home/ritwika/data/1.hc4250_Oct'20/results/systems/200x200_r02_d06.dat", "w");
    fp = fopen("E:\\Stoleriu\\C\\special\\3d\\generare\\2022\\TiOX\\100x100_r02_d06_Cryst10x10.dat", "w");
    for (int i = 0; i < nPart; i++)
    {
        fprintf(fp, "%20.15lf  %20.15lf  %20.15lf  %20.15lf %d\n", Medium[i].x, Medium[i].y, Medium[i].z, Medium[i].r, Medium[i].idxCryst);
    }
    fclose(fp);
    
    

    return 0;
}
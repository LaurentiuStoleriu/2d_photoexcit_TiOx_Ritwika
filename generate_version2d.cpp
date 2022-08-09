//#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>

const int nPart = 1000000;
const int nSide = 1000;

double radius = 0.2;
double distance = 0.6;

struct sReadData
{
    double x, y, z, r;
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
        }
        
    }

    FILE *fp;
    fp = fopen("/home/ritwika/data/1.hc4250_Oct'20/results/systems/200x200_r02_d06.dat", "w");
    //fp = fopen("E:\\Stoleriu\\C\\special\\3d\\generare\\2022\\TiOX\\1000x1000_r02_d06.dat", "w");
    for (int i = 0; i < nPart; i++)
    {
        fprintf(fp, "%20.15lf  %20.15lf  %20.15lf  %20.15lf\n", Medium[i].x, Medium[i].y, Medium[i].z, Medium[i].r);
    }
    fclose(fp);
    
    

    return 0;
}
#include <stdio.h>

const int nPart = 144;
const int nSide = 12;

double radius = 0.2;
double radiusP = 0.15;
double distance = 0.6;

struct sReadData
{
    double x, y, z, r;
};

struct sReadData Medium[nPart];
//struct sReadData Unit[nPartUnit];
int main()
{
   
    double stepX = distance + 2.0 * radius;
    double stepY = distance + 2.0 * radius;
    //double stepP = distance + 2.0 * radius;
    double x0 = radius;
    double y0 = radius;

    for (int i = 1; i < nSide+1; i++)
    {
        for (int j = 1; j < nSide+1; j++)
        {
            Medium[(i-1)*nSide + (j-1)].x = x0 + (j-1) * stepX;
            Medium[(i-1)*nSide + (j-1)].y = y0 + (i-1) * stepY;
            Medium[(i-1)*nSide + (j-1)].z = 0.0;
            //Medium[(i-1)*nSide + (j-1)].r = radius;
            if (i%3 == 0 || j%3 == 0)
                Medium[(i-1)*nSide + (j-1)].r = radiusP;
            else
                Medium[(i-1)*nSide + (j-1)].r = radius;
        }
        
    }

    

    //printf("%20.15lf",Unit[0].x);

    FILE *fp;
    fp = fopen("/home/ritwika/data/1.hc4250_Oct'20/results/systems/12x12_r02_d06_microcrystal.dat", "w");
    for (int i = 0; i < nPart; i++)
    {
        fprintf(fp, "%20.15lf  %20.15lf  %20.15lf  %20.15lf\n", Medium[i].x, Medium[i].y, Medium[i].z, Medium[i].r);
    }
    fclose(fp);
    
    return 0;
}
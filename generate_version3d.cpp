#include <stdio.h>

int main()
{
    const int nPart = 343000;
    const int nSide = 70;

    double radius = 0.2;
    double distance = 0.6;

    struct sReadData
    {
        double x, y, z, r;
    };

    struct sReadData Medium[nPart];

    double stepX = distance + 2.0 * radius;
    double stepY = distance + 2.0 * radius;
    double stepZ = distance + 2.0 * radius;

    double x0 = radius;
    double y0 = radius;
    double z0 = radius;

    for (int i = 0; i < nSide; i++)
    {
        for (int j = 0; j < nSide; j++)
        {
            for (int k = 0; k < nSide; k++)
            {

            
                Medium[i*nSide*nSide + j*nSide + k].x = x0 + k * stepX;
                Medium[i*nSide*nSide + j*nSide + k].y = y0 + j * stepY;
                Medium[i*nSide*nSide + j*nSide + k].z = z0 + i * stepZ;
                Medium[i*nSide*nSide + j*nSide + k].r = radius;
            }
        }
        
    }

    FILE *fp;
    fp = fopen("/home/ritwika/data/1.hc4250_Oct'20/results/systems/70x70x70_r02_d06.dat", "w");
    for (int i = 0; i < nPart; i++)
    {
        fprintf(fp, "%20.15lf  %20.15lf  %20.15lf  %20.15lf\n", Medium[i].x, Medium[i].y, Medium[i].z, Medium[i].r);
    }
    fclose(fp);
    
    

    return 0;
}
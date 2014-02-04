#include "function.h"
#include "stdio.h"
#include "string.h"
#include "portable_time.h"
#include "omp.h"

int main()
{
    double AllTime=PortableGetTime();
    double x0=0.0, y0=0.0, z0=0.0;
    double xn=10.0, yn=10.0, zn=10.0;
    int Sx=300, Sy=300, Sz=300, St=100;
    double * masprev;
    double * masnext;
    masprev=new double[Sx*Sy*Sz];
    masnext=new double[Sx*Sy*Sz];
    double dx=(xn-x0)/Sx, dy=(yn-y0)/Sy, dz=(zn-z0)/Sz;
   
    FILE* filex=fopen("filex.txt","w");
    FILE* filey=fopen("filey.txt","w");
    FILE* filez=fopen("filez.txt","w");

    double dt=0.00001;				//выбираем dt
    omp_set_dynamic(0);      // запретить библиотеке openmp менять число потоков во время исполнения
    omp_set_num_threads(4); // установить число потоков в 10
    memset(masprev, 0, Sx*Sy*Sz*sizeof(double));
	memset(masnext, 0, Sx*Sy*Sz*sizeof(double));
    for (int x=1; x<Sx-1; x++)
        for(int y=1; y<Sy-1; y++)
            for(int z=1; z<Sz-1; z++)
                masprev[x+y*Sx+z*Sx*Sy]=u(x0+dx*x, y0+dy*y, z0+dz*z);
    fprintf(filex,"%e\n", dx);
    fprintf(filey,"%e\n", dy);
    fprintf(filez,"%e\n", dz);
    fprintf(filex,"%i\n", Sx);
    fprintf(filey,"%i\n", Sy);
    fprintf(filez,"%i\n", Sz);
    for(int x=0; x<Sx; x++)
        fprintf(filex,"%lf ", masprev[x+49*Sx+49*Sx*Sy]);
    for(int y=0; y<Sy; y++)
        fprintf(filey,"%lf ", masprev[49+y*Sx+49*Sx*Sy]);
    for(int z=0; z<Sz; z++)
        fprintf(filez,"%lf ", masprev[49+49*Sx+z*Sx*Sy]);
    fprintf(filex,"\n");
    fprintf(filey,"\n");
    fprintf(filez,"\n");

	double Time=PortableGetTime();

    for (int t=1; t<St; t++)
    {        
		#pragma omp parallel for
        for (int z=1; z<Sz-1; z++)
        {
            for(int y=1; y<Sy-1; y++)
            {
                for(int x=1; x<Sx-1; x++)
                {
                    masnext[x+y*Sx+z*Sx*Sy]=
                        dt*((masprev[(x+1)+y*Sx+z*Sx*Sy]-2*masprev[x+y*Sx+z*Sx*Sy]+masprev[(x-1)+y*Sx+z*Sx*Sy])/(dx*dx)
                       +(masprev[x+(y+1)*Sx+z*Sx*Sy]-2*masprev[x+y*Sx+z*Sx*Sy]+masprev[x+(y-1)*Sx+z*Sx*Sy])/(dy*dy)
                       +(masprev[x+y*Sx+(z+1)*Sx*Sy]-2*masprev[x+y*Sx+z*Sx*Sy]+masprev[x+y*Sx+(z-1)*Sx*Sy])/(dz*dz)
                       +f(x0+dx*x, y0+dy*y, z0+dz*z)-masprev[x+y*Sx+z*Sx*Sy])+masprev[x+y*Sx+z*Sx*Sy];
                }
            }
        }

        double* tmp=masprev;
        masprev=masnext;
        masnext=tmp;
    }

    Time=PortableGetTime()-Time;

    fprintf(filex,"%i\n", Sx);
    fprintf(filey,"%i\n", Sy);
    fprintf(filez,"%i\n", Sz);
    for(int x=0; x<Sx; x++)
        fprintf(filex,"%lf ", masprev[x+49*Sx+49*Sx*Sy]);
    for(int y=0; y<Sy; y++)
        fprintf(filey,"%lf ", masprev[49+y*Sx+49*Sx*Sy]);
    for(int z=0; z<Sz; z++)
        fprintf(filez,"%lf ", masprev[49+49*Sx+z*Sx*Sy]);
    fprintf(filex,"\n");
    fprintf(filey,"\n");
    fprintf(filez,"\n");

    AllTime=PortableGetTime()-AllTime;
    printf(" %lf \n %lf \n",Time, AllTime);
    fclose(filex);
    fclose(filey);
    fclose(filez);

	delete[] masprev;
	delete[] masnext;

    return 0;
}
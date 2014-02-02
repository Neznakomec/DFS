#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h> // memcpy
#include <fstream>
#include "mpi.h"
#include "portable_time.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

double t_begin, t_end, time_counter;

double sigma = 1.0;
double k = 1.0;

double TStart, XStart, YStart, ZStart;
double TEnd, XEnd, YEnd, ZEnd;

double** UX;
double** UY;
double** UZ;

double up_f = 4;
//ВРЕМЯ
double GetTime()
{
	return time_counter;
}

double f_begin(double x, double y, double z){ // начальное условие, функция u(0,x,y,z)=ksi(x,y,z)
    double f=0.0;
    if (x>5.0)
        x=0;
    if (y>5.0)
        y=0;
    if (z>5.0)
        z=0;
    f=sin(x*M_PI/5.0)+sin(y*M_PI/5.0)+sin(z*M_PI/5.0);
    return f;
}

double f(double t, double x, double y, double z){ // функция, задающая внешнее воздействие
	double fxyz=0.0;
    if (((x>=5.0) && (x<=5.2)) || ((y>=5.0) && (y<=5.2)) || ((z>=5.0) && (z<=5.2)))
        fxyz=4.0;
    return fxyz;
}

double f_left_x(double time, double y, double z){ // левое граничное условие по x
	return 0;
}

double f_right_x(double time, double y, double z){ // правое граничное условие по x
	return 0;
}

double f_left_y(double time, double x, double z){ // левое граничное условие по y
	return 0;
}

double f_right_y(double time, double x, double z){ // правое граничное условие по y
	return 0;
}
double f_left_z(double time, double x, double y){ // левое граничное условие по z
	return 0;
}

double f_right_z(double time, double x, double y){ // правое граничное условие по z
	return 0;
}

void calculation(	
	int Nx, int Ny, int Nz,	//количество интервалов
	int NTimes     //количество интервалов времени
	)
{
	t_begin = PortableGetTime();
	//
	double XStart = 0.0; double YStart = 0.0; double ZStart = 0.0;
	double TStart = 0.0;
	double deltaX = (XEnd - XStart) / Nx;
	double deltaY = (YEnd - YStart) / Ny;
	double deltaZ = (ZEnd - ZStart) / Nz;
	double deltaT = (TEnd - TStart) / NTimes;

	double ***data[2];
	int stepsX = Nx + 1;
	int stepsY = Ny + 1;
	int stepsZ = Nz + 1;
	int stepsT = NTimes + 1;

	data[0] = new double**[stepsX];
	data[1] = new double**[stepsX];
	for (int i=0; i<stepsX; i++)
	{
		data[0][i] = new double*[stepsY];
		data[1][i] = new double*[stepsY];
		for (int j=0; j<stepsY; j++)
		{
			data[0][i][j] = new double[stepsZ]();
			data[1][i][j] = new double[stepsZ]();
		}
	}

	double curx, cury, curz, curt;

	UX = new double*[stepsT];
	UY = new double*[stepsT];
	UZ = new double*[stepsT];
	for (int i=0; i<stepsT; i++)
	{
		UX[i] = new double[stepsX];
		UY[i] = new double[stepsY];
		UZ[i] = new double[stepsZ];
	}

	int projX, projY, projZ; // моменты времени, по которым мы делаем сечения
	double temp;
	temp = (XStart+XEnd)/2;
	projX =(int)( (temp - XStart) / deltaX );
	temp = (YStart+YEnd)/2;
	projY =(int)( (temp - YStart) / deltaY );
	temp = (ZStart+ZEnd)/2;
	projZ =(int)( (temp - ZStart) / deltaZ );

	// считаем в начальное время:
	for (int i=0; i<stepsX; i++)
		for (int j=0; j<stepsY; j++)
			for (int k=0; k<stepsZ; k++)
			{
				curx = XStart + i*deltaX;
				cury = YStart + j*deltaY;
				curz = ZStart + k*deltaZ;

				data[0][i][j][k] = f_begin(curx, cury, curz);
			};

	for (int i=0; i<stepsX; i++)
		UX[0][i] = data[0][i][projY][projZ];
	for (int j=0; j<stepsY; j++)
		UY[0][j] = data[0][projX][j][projZ];
	for (int k=0; k<stepsZ; k++)
		UZ[0][k] = data[0][projX][projY][k];

	//считаем во все остальные моменты времени
	for(int t = 1; t < stepsT; t++)
	{
		// строим электронные плотности во все моменты времени
		int curr = t % 2;	// вычисление в n+1 момент времени
		int prev = (t-1) % 2; // по n-ному моменту времени

		int i;
		for (i=1; i<stepsX-1; i++)
		{
			for (int j=1; j<stepsY-1; j++)
			{
				for (int k=1; k<stepsZ-1; k++)
				{
					data[curr][i][j][k] = (data[prev][i+1][j][k] - 2*data[prev][i][j][k] + data[prev][i-1][j][k]) / (deltaX*deltaX);
					data[curr][i][j][k] += (data[prev][i][j+1][k] - 2*data[prev][i][j][k] + data[prev][i][j-1][k]) / (deltaY*deltaY);
					data[curr][i][j][k] += (data[prev][i][j][k+1] - 2*data[prev][i][j][k] + data[prev][i][j][k-1]) / (deltaZ*deltaZ);

					data[curr][i][j][k] *= sigma;
					data[curr][i][j][k] -= k*data[prev][i][j][k];
					data[curr][i][j][k] += f(TStart+(t-1)*deltaT, XStart+i*deltaX, YStart+j*deltaY, ZStart+k*deltaZ);
					data[curr][i][j][k] *= deltaT;

					data[curr][i][j][k] += data[prev][i][j][k];
				};
			};
		};


		//считаем граничные условия
		curt = TStart + t*deltaT;

		for (int i=0; i<stepsX; i++)
			for (int j=0; j<stepsY; j++)
			{
				data[curr][i][j][0] = f_left_z(curt, XStart+i*deltaX, YStart+j*deltaY);
				data[curr][i][j][stepsZ-1] = f_right_z(curt, XStart+i*deltaX, YStart+j*deltaY);
			};

		for (int j=0; j<stepsY; j++)
			for (int k=0; k<stepsZ; k++)
			{
				data[curr][0][j][k] = f_left_x(curt, YStart+j*deltaY, ZStart+k*deltaZ);
				data[curr][stepsX-1][j][k] = f_right_x(curt, YStart+j*deltaY, ZStart+k*deltaZ);
			};

		for (int i=0; i<stepsX; i++)
			for (int k=0; k<stepsZ; k++)
			{
				data[curr][i][0][k] = f_left_y(curt, XStart+i*deltaX, ZStart+k*deltaZ);
				data[curr][i][stepsY-1][k] = f_left_y(curt, XStart+i*deltaX, ZStart+k*deltaZ);
			};

		//собираем данные
		for (int i=0; i<stepsX; i++)
			UX[t][i] = data[curr][i][projY][projZ];
		for (int j=0; j<stepsY; j++)
			UY[t][j] = data[curr][projX][j][projZ];
		for (int k=0; k<stepsZ; k++)
			UZ[t][k] = data[curr][projX][projY][k];

	}

	for (int i=0; i<stepsX; i++)
	{
		for (int j=0; j<stepsY; j++)
		{
			delete [](data[0][i][j]);
			delete [](data[1][i][j]);
		};

		delete [](data[0][i]);
		delete [](data[1][i]);
	};
	//
	t_end = PortableGetTime();
	time_counter = t_end - t_begin;
	delete [](data[0]);
	delete [](data[1]);
}

void calc_arguments(int argc, char**argv, int& Nx, int& Ny, int& Nz, int& NTimes)
{
	XStart = YStart = ZStart = TStart = 0.0;

	XEnd = atof(argv[1]); // конец отрезка по x
	YEnd = atof(argv[2]); // конец отрезка по y
	ZEnd = atof(argv[3]); // конец отрезка по z

	Nx = atoi(argv[4]); // количество интервалов, например 100
	Ny = atoi(argv[5]); // количество интервалов, например 100
	Nz = atoi(argv[6]); // количество интервалов, например 100

	double deltaX = (XEnd - XStart) / Nx;
	double deltaY = (YEnd - YStart) / Ny; 
	double deltaZ = (ZEnd - ZStart) / Nz; 

	TStart = 0.0;
	TEnd = atof(argv[7]); // конец времени моделирования
	NTimes = atoi(argv[8]); // количество интервалов по времени, например 100

	double deltaT = (TEnd - TStart) / NTimes;
}

int run(int argc, char **argv)
{
	int Nx, Ny, Nz, NTimes;
	calc_arguments(argc, argv, Nx, Ny, Nz, NTimes);

	int rank, size;
	MPI_Init(NULL, NULL);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	calculation(Nx, Ny, Nz, NTimes);

	MPI_Finalize();
	//вывод в файл
	FILE* filex=fopen("filex.txt","w");
    FILE* filey=fopen("filey.txt","w");
    FILE* filez=fopen("filez.txt","w");

	for (int t = 0; t <= NTimes;  t++)
	{
		fprintf(filex, "time = %i\n", t);
		fprintf(filey, "time = %i\n", t);
		fprintf(filez, "time = %i\n", t);
		for (int i = 0; i <= Nx; i++)
			fprintf(filex, "%lf ", UX[t][i]);
		for (int i = 0; i <= Nx; i++)
			fprintf(filey, "%lf ", UY[t][i]);
	for (int i = 0; i <= Nx; i++)
			fprintf(filez, "%lf ", UZ[t][i]);
        
	fprintf(filex,"\n");
    fprintf(filey,"\n");
    fprintf(filez,"\n");
	}

	fclose(filex);
    fclose(filey);
    fclose(filez);



	return 0;
};

double** getUX()
{
	return UX;
}

double** getUY()
{
	return UY;
}

double** getUZ()
{
	return UZ;
}
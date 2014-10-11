//============================================================================
// Name        : Shallow_Water_OMP.cpp a parallelized OMP version of Shallow_Water_Sequential.cpp
// Authors     : Diego Montufar, Andres Chaves
// Version     : 1.0
// Copyright   : Copyright 2014
// Description : Method-> Finite Difference Method with fourth order central differences
//			     and Fourth Order Runge-Kutta Method
//				 Compile: g++ -fopenmp -o Shallow_Water_Sequential Shallow_Water_Sequential.cpp
//				 Run: ./Shallow_Water_Sequential or just Shallow_Water_Sequential (Windows)
//				 Output-> Assignment1_Files_OMP/Shallow_Water_#.csv numerated files with points on each timestep
//               Visualize-> Open generated csv files with Paraview. 
//============================================================================

#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <limits>
#include <math.h>
#include <omp.h>


using namespace std;

// Global variables
const double	x_min		= 0.00;
const double	x_max		= 100.00;
const double	y_min		= 0.00;
const double	y_max		= 100.00;
const double	t_min		= 0.00;
const double	t_max		= 100.00;
const double	Delta_x		= 1;
const double	Delta_y		= 1;
const double	Delta_t		= 0.1;
const double    g           = 9.81;

const int		N_x			=  (x_max-x_min)/Delta_x+1;
const int		N_y			=  (y_max-y_min)/Delta_y+1;
const int		N_t			=  (t_max-t_min)/Delta_t+1;

// Function declarations
void	f(double** kVx, double** kVy, double** kH, double** phiVx, double** phiVy, double** phiH);
void	write(fstream& file, double** phiH, double** phiVx, double** phiVy,int N_x, int N_y, int t);


int		main(int argc, char** argv)
{
	char myFileName[64];
	fstream	file;

	// Simulation parameters
	double		x			= 0;
	double		y			= 0;
	double		t			= 0;
	int			i			= 0;
	int         j           = 0;
	int			l			= 0;
	double		wtime		= omp_get_wtime();			// Record the starting time
	int			N_Threads	= omp_get_max_threads();

	// Allocate arrays
	cout << "Allocating arrays... " << endl;
		
	double**		tempPhiVx		= new double* [N_x];
	double**		tempPhiVy		= new double* [N_x];
	double**		tempPhiH		= new double* [N_x];

    double**		phiVx			= new double* [N_x];
	double**		phiVy			= new double* [N_x];
	double**		phiH			= new double* [N_x];
	
	double**		k1Vx			= new double* [N_x];
	double**		k1Vy			= new double* [N_x];
	double**		k1H		    	= new double* [N_x];
	double**		k2Vx			= new double* [N_x];
	double**		k2Vy			= new double* [N_x];
	double**		k2H			    = new double* [N_x];
	double**		k3Vx			= new double* [N_x];
	double**		k3Vy			= new double* [N_x];
	double**		k3H			    = new double* [N_x];
	double**		k4Vx			= new double* [N_x];
	double**		k4Vy			= new double* [N_x];
	double**		k4H			    = new double* [N_x];

	for(i=0; i<N_x; i++)
	{
		tempPhiVx[i]	= new double[N_y];
		tempPhiVy[i]	= new double[N_y];
		tempPhiH[i]		= new double[N_y];
		phiVx[i]		= new double[N_y];
		phiVy[i]		= new double[N_y];
		phiH[i]			= new double[N_y];
		k1Vx[i]			= new double[N_y];
		k1Vy[i]			= new double[N_y];
		k1H[i]	    	= new double[N_y];
		k2Vx[i]			= new double[N_y];
		k2Vy[i]			= new double[N_y];
		k2H[i]     		= new double[N_y];
		k3Vx[i]			= new double[N_y];
		k3Vy[i]			= new double[N_y];
		k3H[i]			= new double[N_y];
		k4Vx[i]			= new double[N_y];
		k4Vy[i]			= new double[N_y];
		k4H[i]			= new double[N_y];

	}

	
	/*memset(k1Vx, 0.0, N_x*N_y*sizeof(k1Vx[0][0]));
	memset(k1Vy, 0.0, N_x*N_y*sizeof(k1Vy[0][0]));
	memset(k1H, 0.0, N_x*N_y*sizeof(k1H[0][0]));
	memset(k2Vx, 0.0, N_x*N_y*sizeof(k2Vx[0][0]));
	memset(k2Vy, 0.0, N_x*N_y*sizeof(k2Vy[0][0]));
	memset(k2H, 0.0, N_x*N_y*sizeof(k2H[0][0]));
	memset(k3Vx, 0.0, N_x*N_y*sizeof(k3Vx[0][0]));
	memset(k3Vy, 0.0, N_x*N_y*sizeof(k3Vy[0][0]));
	memset(k3H, 0.0, N_x*N_y*sizeof(k3H[0][0]));
	memset(k4Vx, 0.0, N_x*N_y*sizeof(k4Vx[0][0]));
	memset(k4Vy, 0.0, N_x*N_y*sizeof(k4Vy[0][0]));
	memset(k4H, 0.0, N_x*N_y*sizeof(k4H[0][0]));*/

	// Set initial condition
	cout << "Setting initial conditions... " << endl;
	t	= t_min;
	for(i=0; i<N_x; i++)
	{  
		for(j=0; j<N_y; j++)
		{
			k1Vx[i][j] = k1Vy[i][j] = k1H[i][j] = 0;
			k2Vx[i][j] = k2Vy[i][j] = k2H[i][j] = 0;
			k3Vx[i][j] = k3Vy[i][j] = k3H[i][j] = 0;
			k4Vx[i][j] = k4Vy[i][j] = k4H[i][j] = 0;
			x		= x_min + Delta_x*i;
			y		= y_min + Delta_y*j;

			phiH[i][j]	= 1.00+( 0.5*exp ( (-1.0/25.0)*(pow ( (x-30.0 ),2.0 ) + pow ( ( y-30.0 ),2.0) ) ));
			//cout << fixed << phiH[i][j]<< endl;
			phiVx[i][j] = 0;
			phiVy[i][j] = 0;
		}
	}

	// Write the solution
	sprintf(myFileName, "Assignment1_Files_OMP/Shallow_Water_%d.csv", 1);
	file.open(myFileName, ios::out);
	write(file, phiH, phiVx, phiVy, N_x, N_y, 0);
	file.close();

	cout << "Time Marching Loop.. " << endl;
	
	#pragma omp parallel default(shared) private(i, l, j)
	{
		// Time marching loop
		for(l=0; l<N_t-1; l++)
		{
			#pragma omp single
			{
			   t += Delta_t;
			}
	//		cout << "Calling f(k) " << endl;
			f(k1Vx, k1Vy, k1H, phiVx, phiVy, phiH);
	//		cout << "Called f(k) " << endl;
			#pragma omp for schedule(static)
			for(i=0; i<N_x; i++)
			{
				for(j=0; j<N_y; j++)
				{
					tempPhiVx[i][j]	= phiVx[i][j] + (Delta_t/2)*k1Vx[i][j];
					tempPhiVy[i][j]	= phiVy[i][j] + (Delta_t/2)*k1Vy[i][j];
					tempPhiH[i][j]	= phiH[i][j] + (Delta_t/2)*k1H[i][j];
				}
			}
			//cout << "Calling f(k2) " << endl;
		    f(k2Vx, k2Vy, k2H, tempPhiVx, tempPhiVy, tempPhiH);

			#pragma omp for schedule(static)
			for(i=0; i<N_x; i++)
			{
				for(j=0; j<N_y; j++)
				{
					tempPhiVx[i][j]	= phiVx[i][j] + (Delta_t/2)*k2Vx[i][j];
					tempPhiVy[i][j]	= phiVy[i][j] + (Delta_t/2)*k2Vy[i][j];
					tempPhiH[i][j]	= phiH[i][j] + (Delta_t/2)*k2H[i][j];
				}
			}
		
			f(k3Vx, k3Vy, k3H, tempPhiVx, tempPhiVy, tempPhiH);

			#pragma omp for schedule(static)
			for(i=0; i<N_x; i++)
			{
				for(j=0; j<N_y; j++)
				{
					tempPhiVx[i][j]	= phiVx[i][j] + Delta_t*k3Vx[i][j];
					tempPhiVy[i][j]	= phiVy[i][j] + Delta_t*k3Vy[i][j];
					tempPhiH[i][j]	= phiH[i][j] + Delta_t*k3H[i][j];
				}			
			}

			f(k4Vx, k4Vy, k4H, tempPhiVx, tempPhiVy, tempPhiH);

			#pragma omp for schedule(static)
			for(i=0; i<N_x; i++)
			{
				for(j=0; j<N_y; j++)
				{
					phiVx[i][j]	=  phiVx[i][j] + Delta_t*(k1Vx[i][j]/6 + k2Vx[i][j]/3 + k3Vx[i][j]/3 + k4Vx[i][j]/6);
					phiVy[i][j]	=  phiVy[i][j] + Delta_t*(k1Vy[i][j]/6 + k2Vy[i][j]/3 + k3Vy[i][j]/3 + k4Vy[i][j]/6);
					phiH[i][j]	=  phiH[i][j] + Delta_t*(k1H[i][j]/6 + k2H[i][j]/3 + k3H[i][j]/3 + k4H[i][j]/6);
					//if (phiH[i][j] != phiH[i][j])
						//cout << phiH[i][j] << "," << k1H[i][j] << "," << k2H[i][j] << "," << k3H[i][j] << "," << k4H[i][j]<<endl;
				}
			}

			#pragma omp single
			{
            	// Write the solution
				sprintf(myFileName, "Assignment1_Files_OMP/Shallow_Water_%d.csv", l+1);
				file.open(myFileName, ios::out);
				write(file, phiH, phiVx, phiVy, N_x, N_y, t+1);
				file.close();
			}
		
			//Print Time marching loop
			cout << "t = " << t << endl;
		}
	}
	
	file.close();

	wtime	= omp_get_wtime() - wtime;	// Record the end time and calculate elapsed time
	cout << "Simulation took " << wtime/N_t << " seconds per time step with " << N_Threads << " threads" << endl;

	// Deallocate arrays
	delete [] tempPhiVx;
	delete [] tempPhiVy;
	delete [] tempPhiH;
	delete [] k1Vx;
	delete [] k1Vy;
	delete [] k1H;	
	delete [] k2Vx;
	delete [] k2Vy;
	delete [] k2H;	
	delete [] k3Vx;
	delete [] k3Vy;
	delete [] k3H;	
	delete [] k4Vx;
	delete [] k4Vy;
	delete [] k4H;	

	delete [] phiVx;
	delete [] phiVy;
	delete [] phiH;
	
	return 0;
}


void	f(double** kVx, double** kVy, double** kH, double** phiVx, double** phiVy, double** phiH)
{
	//cout << "Function f " << endl;
	int ip1, ip2, ip3, ip4, jp1, jp2, jp3, jp4;
	#pragma omp for
	for(int i=0; i<N_x; i++)
	{
		for(int j=0; j<N_y; j++)
		{
			//cout << "Function f " << i << "," << j << endl;
			
			ip1=i+1;
			ip2=i-1;
			ip3=i+2;
			ip4=i-2;
			jp1=j+1;
			jp2=j-1;
			jp3=j+2;
			jp4=j-2;
            
            if(i==N_x-1){
                 ip1=0;
                 ip3=1;
			}
            
            if(i==N_x-2){
                 ip3=0;
			}
            
            if(i==0){
                ip2=N_x-1;
                ip4=N_x-2;
			}
            
            if(i==1){
                ip4=N_x-1;
			}
            
            if(j==N_y-1){
                  jp1=0;
                  jp3=1;
			}
            
            if(j==N_y-2){
                 jp3=0;
			}
            
            if(j==0){
                jp2=N_y-1;
                jp4=N_y-2;
			}
            
            if(j==1){
                jp4=N_y-1;
			}
			//cout << "Border condition assesed, calculating Ks " << phiH[ip2][j] << endl;
            kVx[i][j]= (-g/(12*Delta_x)*(phiH[ip4][j]-8*(phiH[ip2][j])+8*phiH[ip1][j]-phiH[ip3][j])) - (phiVx[i][j]/(12*Delta_x)*(phiVx[ip4][j]-8*(phiVx[ip2][j])+8*phiVx[ip1][j]-phiVx[ip3][j])) - (phiVx[i][j]/(12*Delta_y)*(phiVx[i][jp4]-8*(phiVx[i][jp2])+8*phiVx[i][jp1]-phiVx[i][jp3]));
            kVy[i][j]= (-g/(12*Delta_y)*(phiH[i][jp4]-8*(phiH[i][jp2])+8*phiH[i][jp1]-phiH[i][jp3])) - (phiVx[i][j]/(12*Delta_x)*(phiVy[ip4][j]-8*(phiVy[ip2][j])+8*phiVy[ip1][j]-phiVy[ip3][j])) - (phiVy[i][j]/(12*Delta_y)*(phiVy[i][jp4]-8*(phiVy[i][jp2])+8*phiVy[i][jp1]-phiVy[i][jp3]));
            kH[i][j] = (-phiVx[i][j]/(12*Delta_x)*(phiH[ip4][j]-8*(phiH[ip2][j])+8*phiH[ip1][j]-phiH[ip3][j])) - (phiH[i][j]/(12*Delta_x)*(phiVx[ip4][j]-8*(phiVx[ip2][j])+8*phiVx[ip1][j]-phiVx[ip3][j])) - (phiVy[i][j]/(12*Delta_y)*(phiH[i][jp4]-8*(phiH[i][jp2])+8*phiH[i][jp1]-phiH[i][jp3])) - (phiH[i][j]/(12*Delta_y)*(phiVy[i][jp4]-8*(phiVy[i][jp2])+8*phiVy[i][jp1]-phiVy[i][jp3]));

	//		if (kVx[i][j] == 0 || kVy[i][j] == 0 || kH[i][j] == 0)
	  // 			cout << kVx[i][j] << "," << kVy[i][j] << "," << kH[i][j] << endl;

		}
	}
	return;
}

void	write(fstream& file, double** phiH, double** phiVx, double** phiVy, int N_x, int N_y, int t)
{
	for(int i=0; i<N_x; i++)
	{
		double coordX = i*Delta_x;
		for(int j=0; j<N_y; j++)
		{
			double coordY = j*Delta_y;
			//file << std::fixed << std::setprecision(15) << phiH[i][j] << "\t"; //One single file
			file << std::fixed << std::setprecision(4) << phiH[i][j] <<"," << coordX << "," << coordY<< "\n"; //multiple files
		}
		//file << endl;
	}
	return;
}

//============================================================================
// Name        : Shallow_Water_MPI.cpp a parallelized MPI version of Shallow_Water_Sequential.cpp 
// Authors     : Diego Montufar, Andres Chaves
// Version     : 1.0
// Copyright   : Copyright 2014
// Description : Method-> Finite Difference Method with fourth order central differences
//			     and Fourth Order Runge-Kutta Method
//				 Compile: mpicxx -o Shallow_Water_MPI Shallow_Water_MPI.cpp
//				 Run: mpirun -np 4 ./Shallow_Water_MPI
//				 Output-> Assignment1_Files_MPI/Shallow_Water_MPI_#.csv numerated files with points on each timestep
//               Visualize-> Open generated csv files with Paraview. 
//============================================================================
#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <limits>
#include <mpi.h>
#include <math.h>

using namespace std;

const int 			N_D = 2;
const int 			X = 0;
const int 			Y = 1;
const double    	t_min  =  0.00;
const double    	t_max  = 100.00;
const int 			numElementsPerBlock	= 1;
const double    	Delta_t =  0.1;
const int 			N_t = (t_max-t_min)/Delta_t;
const double    	g = 9.81;


// Function declarations
void	f(double** kVx, double** kVy, double** kH, double** phiVx, double** phiVy, double** phiH, int myN_x, int myN_y, double delta_x, double delta_y, int flag);
void	write(fstream& file, double** phiH, double** phiVx, double** phiVy,int N_x, int N_y, int t, int* myCoords,double delta_x, double delta_y);
void	exchange(double** phiH, double** phiVx, double** phiVy , int myN_x, int myN_y, int myID,
	             int rightNeihgbor, int leftNeighbor, int topNeighbor, int bottomNeighbor, MPI_Comm comm2D, MPI_Datatype strideType);

int	main(int argc, char** argv){
	
	int				myID;
	int				N_Procs;
  
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &N_Procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myID);

	int				N				= sqrt(N_Procs);
	int				dimensions[N_D]	= {N, N};	
	int				isPeriodic[N_D]	= {1, 1};
	int				myCoords[N_D]	= {0, 0};
	int				myN_x			= 200;
	int				myN_y			= 200;
	int				N_x				= myN_x*dimensions[X];
	int				N_y				= myN_y*dimensions[Y];
	int				myiStart		= 0;
	int				myiEnd			= 0;
	int				myjStart		= 0;
	int				myjEnd			= 0;
	double			x_min			= 0.00;
	double			x_max			= 100.00;
	double			y_min			= 0.00;
	double			y_max			= 100.00;
	double			Delta_x			= ((double)x_max-(double)x_min)/(double)(N_x);
	double			Delta_y			= ((double)y_max-(double)y_min)/(double)(N_y);
	double			Delta_xy		= Delta_x;
	int				i				= 0;
	int				j				= 0;
	int				l				= 0;
	int				leftNeighbor	= 0;
	int				rightNeighbor	= 0;
	int				bottomNeighbor	= 0;
	int				topNeighbor		= 0;
	int				reorder			= 1;
	double			wtime			= 0;
	double          x               = 0;
    double          y               = 0;
    double          t               = 0;
	fstream			file;
	char			myFileName[64];
	MPI_Status		status;
	MPI_Datatype	strideType;
	MPI_Comm		Comm2D;

	
	if(myID==0)	{
		wtime	= MPI_Wtime();
	}
	
	if ( myID == 0 ) cout << "There are " << N_Procs << " procs. Distributing " << myN_x << " x " << myN_y << " elements of"<<N_x <<" x" << N_y << "total. Delta x:" << Delta_x <<" Delta y:"<< Delta_y  << endl;

	// Allocate arrays
	cout << "Allocating arrays.. " << endl;

	// Allocate arrays	
	double** phi				= new double*	[myN_x+4];
	double** r					= new double*	[myN_x+4];
	
    phi[0]						= new double	[(myN_x+4)*(myN_y+4)];
    r[0]						= new double	[(myN_x+4)*(myN_y+4)];
	
	//Array of pointers
	double**		tempPhiVx		= new double* [myN_x+4];
	double**		tempPhiVy		= new double* [myN_x+4];
	double**		tempPhiH		= new double* [myN_x+4];
    double**		phiVx			= new double* [myN_x+4]; //2d allocated array +4 rows for the halo
	double**		phiVy			= new double* [myN_x+4]; //2d allocated array +4 rows for the halo
	double**		phiH			= new double* [myN_x+4]; //2d allocated array +4 rows for the halo
	double**		k1Vx			= new double* [myN_x+4];
	double**		k1Vy			= new double* [myN_x+4];
	double**		k1H		    	= new double* [myN_x+4];
	double**		k2Vx			= new double* [myN_x+4];
	double**		k2Vy			= new double* [myN_x+4];
	double**		k2H			    = new double* [myN_x+4];
	double**		k3Vx			= new double* [myN_x+4];
	double**		k3Vy			= new double* [myN_x+4];
	double**		k3H			    = new double* [myN_x+4];
	double**		k4Vx			= new double* [myN_x+4];
	double**		k4Vy			= new double* [myN_x+4];
	double**		k4H			    = new double* [myN_x+4];

	//2d Array in 1d larger consecutive array
	tempPhiVx[0] = new double [(myN_x+4)*(myN_y+4)];
	tempPhiVy[0] = new double [(myN_x+4)*(myN_y+4)];
	tempPhiH[0] = new double [(myN_x+4)*(myN_y+4)];
	phiVx[0] = new double [(myN_x+4)*(myN_y+4)];
	phiVy[0] = new double [(myN_x+4)*(myN_y+4)];
	phiH[0] = new double [(myN_x+4)*(myN_y+4)];
	k1Vx[0] = new double [(myN_x+4)*(myN_y+4)];
	k1Vy[0] = new double [(myN_x+4)*(myN_y+4)];
	k1H[0] = new double [(myN_x+4)*(myN_y+4)];
	k2Vx[0] = new double [(myN_x+4)*(myN_y+4)];
	k2Vy[0] = new double [(myN_x+4)*(myN_y+4)];
	k2H[0] = new double [(myN_x+4)*(myN_y+4)];
	k3Vx[0] = new double [(myN_x+4)*(myN_y+4)];
	k3Vy[0] = new double [(myN_x+4)*(myN_y+4)];
	k3H[0] = new double [(myN_x+4)*(myN_y+4)];
	k4Vx[0] = new double [(myN_x+4)*(myN_y+4)];
	k4Vy[0] = new double [(myN_x+4)*(myN_y+4)];
	k4H[0] = new double [(myN_x+4)*(myN_y+4)];

	//Divide the array in a 2D way
    for(int i=1, ii=myN_y+4; i<myN_x+4; i++, ii+=(myN_y+4))
	{
		tempPhiVx[i] = &tempPhiVx[0][ii];
		tempPhiVy[i] = &tempPhiVy[0][ii];
		tempPhiH[i]	 = &tempPhiH[0][ii];

		phiVx[i] = &phiVx[0][ii];
		phiVy[i] = &phiVy[0][ii];
		phiH[i]	 = &phiH[0][ii];

		k1Vx[i] = &k1Vx[0][ii];
		k1Vy[i] = &k1Vy[0][ii];
		k1H[i]	 = &k1H[0][ii];

		k2Vx[i] = &k2Vx[0][ii];
		k2Vy[i] = &k2Vy[0][ii];
		k2H[i]	 = &k2H[0][ii];

		k3Vx[i] = &k3Vx[0][ii];
		k3Vy[i] = &k3Vy[0][ii];
		k3H[i]	 = &k3H[0][ii];

		k4Vx[i] = &k4Vx[0][ii];
		k4Vy[i] = &k4Vy[0][ii];
		k4H[i]	 = &k4H[0][ii];

	}
	
	// Map the processors to the N x N grid, determine the new rank and determine the neighboring processes
	MPI_Cart_create(MPI_COMM_WORLD, N_D, dimensions, isPeriodic, reorder, &Comm2D);
	MPI_Comm_rank(Comm2D, &myID);
	MPI_Cart_coords(Comm2D, myID, N_D, myCoords);
	MPI_Cart_shift(Comm2D, X, 1, &leftNeighbor,   &rightNeighbor);
	MPI_Cart_shift(Comm2D, Y, 1, &bottomNeighbor, &topNeighbor);

	myiStart	= myCoords[X]==0   ? 2     : 1;
	myiEnd		= myCoords[X]==N-1 ? myN_x : myN_x+1;
	myjStart	= myCoords[Y]==0   ? 2     : 1;
	myjEnd		= myCoords[Y]==N-1 ? myN_y : myN_y+1;

	//cout << "Process " << myID << " coords (" << myCoords[X] << "," << myCoords[Y] << ")" << "iStart = " << myiStart <<  ", iEnd = " << myiEnd << ", jStart = " << myjStart <<  ", jEnd = " << myjEnd << "X: "<<myN_x <<" " << "Y:"<<myN_y << endl;
	//cout << "Process " << myID << " coords (" << myCoords[X] << "," << myCoords[Y] << ")" << "X: "<<myN_x <<" " << "Y:"<<myN_y << " LN:" << leftNeighbor << " RN:" << rightNeighbor << " TN:" << topNeighbor << " BN:" << bottomNeighbor << endl;

	MPI_Barrier(Comm2D);
	
	// Set initial condition
	cout << myID <<": Setting initial conditions " << endl;
	t = t_min;
	
	for(i=2; i<=myN_x+1; i++)
	{  
	    for(j=2; j<=myN_y+1; j++)
		{
//			cout << myID << "Coord "<< i << "," << j << endl;
			int ii	= myCoords[X]*myN_x + i - 2;
			int jj	= myCoords[Y]*myN_y + j - 2;

			k1Vx[i][j] = k1Vy[i][j] = k1H[i][j] = 0;
			k2Vx[i][j] = k2Vy[i][j] = k2H[i][j] = 0;
			k3Vx[i][j] = k3Vy[i][j] = k3H[i][j] = 0;
			k4Vx[i][j] = k4Vy[i][j] = k4H[i][j] = 0;
			x	= x_min + Delta_x*ii;
			y	= y_min + Delta_y*jj;

			phiH[i][j]	= 1.00+( 0.5*exp ( (-1.0/25.0)*(pow ( (x-30.0 ),2.0 ) + pow ( ( y-30.0 ),2.0) ) ));
			phiVx[i][j] = 0;
			phiVy[i][j] = 0;
		}
	}
	cout << myID <<": Conditions sat " << endl;

	// Create a new datatype to store values on an x boundary
	MPI_Type_vector(myN_x, numElementsPerBlock, myN_y+4, MPI_DOUBLE, &strideType);
	MPI_Type_commit(&strideType);

    // Write the solution
	sprintf(myFileName, "Assignment1_Files_MPI/Shallow_Water_MPI.csv.%d", myID);
	file.open(myFileName, ios::out);
    write(file, phiH, phiVx, phiVy, myN_x, myN_y, 0, myCoords,Delta_x, Delta_y);
	
	cout << myID <<": Time Marching Loop " << endl;
	
	// Time marching loop
	for(l=0; l<N_t-1; l++)
	{
    	t += Delta_t;
		
		// Exchange data  
		exchange(phiH, phiVx, phiVy , myN_x, myN_y, myID, rightNeighbor, leftNeighbor, topNeighbor,  bottomNeighbor, Comm2D, strideType);
		f(k1Vx, k1Vy, k1H, phiVx, phiVy, phiH, myN_x, myN_y,Delta_x, Delta_y,1);

		for(i=2; i<=myN_x+1; i++)
		{
			for(j=2; j<=myN_y+1; j++)
			{
				tempPhiVx[i][j]	= phiVx[i][j] + (Delta_t/2)*k1Vx[i][j];
				tempPhiVy[i][j]	= phiVy[i][j] + (Delta_t/2)*k1Vy[i][j];
				tempPhiH[i][j]	= phiH[i][j] + (Delta_t/2)*k1H[i][j];
			}
		}

		// Exchange data   
		exchange(tempPhiH, tempPhiVx, tempPhiVy , myN_x, myN_y, myID, rightNeighbor, leftNeighbor, topNeighbor,  bottomNeighbor, Comm2D, strideType);
        f(k2Vx, k2Vy, k2H, tempPhiVx, tempPhiVy, tempPhiH, myN_x, myN_y,Delta_x, Delta_y,2);
		
		for(i=2; i<=myN_x+1; i++)
		{
			for(j=2; j<=myN_y+1; j++)
			{
				tempPhiVx[i][j]	= phiVx[i][j] + (Delta_t/2)*k2Vx[i][j];
				tempPhiVy[i][j]	= phiVy[i][j] + (Delta_t/2)*k2Vy[i][j];
				tempPhiH[i][j]	= phiH[i][j] + (Delta_t/2)*k2H[i][j];
			}
		}

		
		// Exchange data  
		exchange(tempPhiH, tempPhiVx, tempPhiVy , myN_x, myN_y, myID, rightNeighbor, leftNeighbor, topNeighbor,  bottomNeighbor, Comm2D, strideType);
		f(k3Vx, k3Vy, k3H, tempPhiVx, tempPhiVy, tempPhiH, myN_x, myN_y,Delta_x, Delta_y,3);
		
		for(i=2; i<=myN_x+1; i++)
		{
			for(j=2; j<=myN_y+1; j++)
			{
				tempPhiVx[i][j]	= phiVx[i][j] + Delta_t*k3Vx[i][j];
				tempPhiVy[i][j]	= phiVy[i][j] + Delta_t*k3Vy[i][j];
				tempPhiH[i][j]	= phiH[i][j] + Delta_t*k3H[i][j];
				if (l==0 && i==14 && j==58)
					cout << fixed << std::setprecision(10)<< myID <<"K333:"<<k3H[i][j]<<endl;
			}			
		}

		// Exchange data  
		exchange(tempPhiH, tempPhiVx, tempPhiVy , myN_x, myN_y, myID, rightNeighbor, leftNeighbor, topNeighbor,  bottomNeighbor, Comm2D, strideType);
		f(k4Vx, k4Vy, k4H, tempPhiVx, tempPhiVy, tempPhiH, myN_x, myN_y,Delta_x, Delta_y,4 );

		for(i=2; i<=myN_x+1; i++)
		{
			for(j=2; j<=myN_y+1; j++)
			{
				phiVx[i][j]	=  phiVx[i][j] + Delta_t*(k1Vx[i][j]/6 + k2Vx[i][j]/3 + k3Vx[i][j]/3 + k4Vx[i][j]/6);
				phiVy[i][j]	=  phiVy[i][j] + Delta_t*(k1Vy[i][j]/6 + k2Vy[i][j]/3 + k3Vy[i][j]/3 + k4Vy[i][j]/6);
				phiH[i][j]	=  phiH[i][j] + Delta_t*(k1H[i][j]/6 + k2H[i][j]/3 + k3H[i][j]/3 + k4H[i][j]/6);
			}
		}
		
		// Write the solution
	    write(file, phiH, phiVx, phiVy, myN_x, myN_y, l+1, myCoords,Delta_x, Delta_y);
    		
		if (myID==0) cout << "t = " << t << endl;
	}
	file.close();
	
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

	if(myID==0)
	{
		wtime	= MPI_Wtime() - wtime;	// Record the end time and calculate elapsed time
		cout << "Simulation took " << wtime/N_t << " seconds per time step with " << N_Procs << " processes" << endl;
	}

    MPI_Finalize();
	return 0;
}


void	f(double** kVx, double** kVy, double** kH, double** phiVx, double** phiVy, double** phiH, int myN_x, int myN_y, double Delta_x, double Delta_y, int flag)
{
	//cout << "Function f " << endl;
	int ip1, ip2, ip3, ip4, jp1, jp2, jp3, jp4;
	
	for(int i=2; i<=myN_x+1; i++)
	{
		for(int j=2; j<=myN_y+1; j++)
		{
			//cout << "Function f " << i << "," << j << endl;
			
			ip1=i+1;ip2=i-1;ip3=i+2;ip4=i-2;jp1=j+1;jp2=j-1;jp3=j+2;jp4=j-2;
            
            kVx[i][j]= (-g/(12*Delta_x)*(phiH[ip4][j]-8*(phiH[ip2][j])+8*phiH[ip1][j]-phiH[ip3][j])) - (phiVx[i][j]/(12*Delta_x)*(phiVx[ip4][j]-8*(phiVx[ip2][j])+8*phiVx[ip1][j]-phiVx[ip3][j])) - (phiVx[i][j]/(12*Delta_y)*(phiVx[i][jp4]-8*(phiVx[i][jp2])+8*phiVx[i][jp1]-phiVx[i][jp3]));
            kVy[i][j]= (-g/(12*Delta_y)*(phiH[i][jp4]-8*(phiH[i][jp2])+8*phiH[i][jp1]-phiH[i][jp3])) - (phiVx[i][j]/(12*Delta_x)*(phiVy[ip4][j]-8*(phiVy[ip2][j])+8*phiVy[ip1][j]-phiVy[ip3][j])) - (phiVy[i][j]/(12*Delta_y)*(phiVy[i][jp4]-8*(phiVy[i][jp2])+8*phiVy[i][jp1]-phiVy[i][jp3]));
            kH[i][j] = (-phiVx[i][j]/(12*Delta_x)*(phiH[ip4][j]-8*(phiH[ip2][j])+8*phiH[ip1][j]-phiH[ip3][j])) - (phiH[i][j]/(12*Delta_x)*(phiVx[ip4][j]-8*(phiVx[ip2][j])+8*phiVx[ip1][j]-phiVx[ip3][j])) - (phiVy[i][j]/(12*Delta_y)*(phiH[i][jp4]-8*(phiH[i][jp2])+8*phiH[i][jp1]-phiH[i][jp3])) - (phiH[i][j]/(12*Delta_y)*(phiVy[i][jp4]-8*(phiVy[i][jp2])+8*phiVy[i][jp1]-phiVy[i][jp3]));
		}
	}
	return;
}

void exchange(double** phiH, double** phiVx, double** phiVy , int myN_x, int myN_y, int myID, 
              int rightNeighbor, int leftNeighbor, int topNeighbor, int bottomNeighbor, MPI_Comm Comm2D, MPI_Datatype strideType){
	MPI_Status		status;

	MPI_Sendrecv(&(phiH[2][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0,	&(phiH[myN_x+2][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiH[3][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0,	&(phiH[myN_x+3][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiH[myN_x][2]),		myN_y,	MPI_DOUBLE,	rightNeighbor,		0,	&(phiH[0][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiH[myN_x+1][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0,	&(phiH[1][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiH[2][3]),			1,		strideType, bottomNeighbor,		0,	&(phiH[2][myN_x+3]),	1,		strideType,	topNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiH[2][2]),			1,		strideType, bottomNeighbor,		0,	&(phiH[2][myN_x+2]),	1,		strideType,	topNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiH[2][myN_x]),		1,		strideType, topNeighbor,		0,	&(phiH[2][0]),			1,		strideType,	bottomNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiH[2][myN_x+1]),	1,		strideType, topNeighbor,		0,	&(phiH[2][1]),			1,		strideType,	bottomNeighbor,		0, Comm2D, &status);

	MPI_Sendrecv(&(phiVx[2][2]),		myN_y,	MPI_DOUBLE,	leftNeighbor,		0,	&(phiVx[myN_x+2][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVx[3][2]),		myN_y,	MPI_DOUBLE,	leftNeighbor,		0,	&(phiVx[myN_x+3][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVx[myN_x][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0,	&(phiVx[0][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVx[myN_x+1][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0,	&(phiVx[1][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVx[2][3]),		1,		strideType, bottomNeighbor,		0,	&(phiVx[2][myN_x+3]),	1,		strideType,	topNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVx[2][2]),		1,		strideType, bottomNeighbor,		0,	&(phiVx[2][myN_x+2]),	1,		strideType,	topNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVx[2][myN_x]),	1,		strideType, topNeighbor,		0,	&(phiVx[2][0]),			1,		strideType,	bottomNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVx[2][myN_x+1]),	1,		strideType, topNeighbor,		0,	&(phiVx[2][1]),			1,		strideType,	bottomNeighbor,		0, Comm2D, &status);

	MPI_Sendrecv(&(phiVy[2][2]),		myN_y,	MPI_DOUBLE,	leftNeighbor,		0,	&(phiVy[myN_x+2][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVy[3][2]),		myN_y,	MPI_DOUBLE,	leftNeighbor,		0,	&(phiVy[myN_x+3][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVy[myN_x][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0,	&(phiVy[0][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVy[myN_x+1][2]),	myN_y,	MPI_DOUBLE,	rightNeighbor,		0,	&(phiVy[1][2]),			myN_y,	MPI_DOUBLE,	leftNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVy[2][3]),		1,		strideType, bottomNeighbor,		0,	&(phiVy[2][myN_x+3]),	1,		strideType,	topNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVy[2][2]),		1,		strideType, bottomNeighbor,		0,	&(phiVy[2][myN_x+2]),	1,		strideType,	topNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVy[2][myN_x]),	1,		strideType, topNeighbor,		0,	&(phiVy[2][0]),			1,		strideType,	bottomNeighbor,		0, Comm2D, &status);
	MPI_Sendrecv(&(phiVy[2][myN_x+1]),	1,		strideType, topNeighbor,		0,	&(phiVy[2][1]),			1,		strideType,	bottomNeighbor,		0, Comm2D, &status);
	
}	

void write(fstream& file, double** phiH, double** phiVx, double** phiVy, int N_x, int N_y, int t, int* myCoords, double Delta_x, double Delta_y)
{
	for(int i=2; i<=N_x+1; i++)
	{ 
	    double coordX = (myCoords[X]*N_x + (i-2))*Delta_x;
		for(int j=2; j<=N_y+1; j++)
		{
			double coordY = (myCoords[Y]*N_y + (j-2))*Delta_y;
			if (t==1 && coordX==6 && coordY==28 )
				//cout << fixed << std::setprecision(10) << coordX << "," << coordY << "=" << phiH[i][j] <<"-"<<i<<","<<j<< endl;
				file << t <<","<< std::fixed << std::setprecision(10) << coordX << "," << coordY<< "," << phiH[i][j] <<"\n";
		}
	}
	return;
}

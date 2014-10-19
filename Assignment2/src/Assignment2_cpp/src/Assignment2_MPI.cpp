/////////////////////////////////////////////////////////////////////////////
//
// Applied Numerical Methods
//
// Assignment 2: The Heat Equation
//
// Group members: Diego Montufar/Andres Chaves
//
// Problem:	rho C dphi/dt = -k Grad^2 phi
//
// Method:		Finite Element Method with linear 2D triangular elements
//				and Implicit Euler Method and Conjugate Gradient Method
//
// Compilation:	mpicxx Assignment2_MPI.cpp -o Assignment2_MPI.exe
//
// Execution:	./Assignment2_MPI /grids/box/Box.grid
//
/////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <mpi.h>
#include "SparseMatrix.h"
#include "Boundary.h"

using namespace std;

// Global variables
const double t_min 	= 0.00;
const double t_max 	= 100.00;
const double Delta_t 	= 0.01;

const double rho 		= 8954.00;
const double C 		= 380.00;
const double k 		= 386.00;
const double h 		= 100.00;
const double Tair 		= 300.00;
const double Qcpu 		= 40000.00;

const int N_t 			= static_cast<int>((t_max - t_min) / Delta_t + 1);
double* buffer 			= NULL; //Added
int bufferSize 			= 0; //Added

// Function declarations
void	readData(char* filename, double**& Points, int**& Faces, int**& Elements, Boundary*& Boundaries, int& myN_p, int& myN_f, int& myN_e, int& myN_b, bool*& yourPoints, int myID);
void	writeData(double* phi, double**& Points, int**& Elements, int& myN_p, int& myN_e, int l, int myID);
void	exchangeData(double* v, Boundary* Boundaries, int myN_b);
void	assembleSystem(SparseMatrix& M, SparseMatrix& K, double* s, double* phi, bool* Free, bool* Fixed, double** Points, int** Faces, int** Elements, Boundary* Boundaries, int myN_p, int myN_f, int myN_e, int myN_b, int myID);
void	solve(SparseMatrix& A, double* phi, double* b, bool* Free, bool* Fixed, Boundary* Boundaries, bool* yourPoints, int myN_b, int myID);
double	computeInnerProduct(double* v1, double* v2, bool* Free, bool* yourPoints, int N_row);

int     main(int argc, char** argv)
{
    // Simulation parameters
    double**		Points		= NULL;
    int**			Faces		= NULL;
    int**			Elements	= NULL;
    Boundary*		Boundaries	= NULL;
    bool*			yourPoints	= NULL; //Added
    double*			buffer		= NULL; //Added
    int				myN_p		= 0; //Changed
    int				myN_f		= 0; //Changed
    int				myN_e		= 0; //Changed
    int				myN_b		= 0; //Changed
    int				myID		= 0; //Added
    int				N_Procs		= 0; //Added
	double			wtime; //Added

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,	&myID);
	MPI_Comm_size(MPI_COMM_WORLD,	&N_Procs);

	if(argc<2)
	{
		if(myID==0)
		{
			cerr << "No grid file specified" << endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
	{
		readData(argv[1], Points, Faces, Elements, Boundaries, myN_p, myN_f, myN_e, myN_b, yourPoints, myID);
	}

	// Allocate arrays
    double*			phi			= new double[myN_p];
    double*			s			= new double[myN_p];
    double*			b			= new double[myN_p];
    bool*			Free        = new bool	[myN_p];
    bool*			Fixed       = new bool	[myN_p];
    double*			AphiFixed	= new double[myN_p];
    SparseMatrix	M;
    SparseMatrix	K;
    SparseMatrix	A;

    if(myID==0)
	{
		wtime	= MPI_Wtime(); //Added
	}

    // Set initial condition
    for(int m=0; m<myN_p; m++)
    {
        phi[m]	= 300;
    }

    assembleSystem(M, K, s, phi, Free, Fixed, Points, Faces, Elements, Boundaries, myN_p, myN_f, myN_e, myN_b, myID); //Changed

    A = M;
    A.subtract(Delta_t, K); // At this point we have A = M-Delta_t*K

    // Compute the column vector to subtract from the right hand side to take account of fixed nodes
    A.multiply(AphiFixed, phi, Free, Fixed);
    exchangeData(AphiFixed, Boundaries, myN_b); //Added
    exchangeData(s, Boundaries, myN_b); //Added

    writeData(phi, Points, Elements, myN_p, myN_e, 0, myID);

    // Time marching loop
    for(int l=0; l<N_t-1; l++)
    {
    	if(myID==0)
    	{
			cout << "t = " << l*Delta_t;
		}

        // Assemble b
        M.multiply(b, phi); // b = M*phi^l
        exchangeData(b, Boundaries, myN_b);
        for(int m=0; m<myN_p; m++)
        {
            b[m]	+= Delta_t*s[m] - AphiFixed[m];
        } // b = M*phi^l + Delta_t*s - A_free,fixed*phi_fixed

        // Solve the linear system
        solve(A, phi, b, Free, Fixed, Boundaries, yourPoints, myN_b, myID);//Changed

        // Write the solution
        if(l%10==0)
		{
			writeData(phi, Points, Elements, myN_p, myN_e, l, myID);
		}
    }

    if(myID==0)
	{
		wtime	= MPI_Wtime() - wtime;	// Record the end time and calculate elapsed time
		cout << "Simulation took " << wtime << " seconds with " << N_Procs << " processes" << endl;
	} //Added

    MPI_Buffer_detach(&buffer, &bufferSize); //Added

    // Deallocate arrays
    for(int boundary=0; boundary<myN_b; boundary++)
    {
        delete [] Boundaries[boundary].indices_;
    }
    delete [] Points[0];
    delete [] Points;
    delete [] Faces[0];
    delete [] Faces;
    delete [] Elements[0];
    delete [] Elements;
    delete [] Boundaries;
    delete [] phi;
    delete [] s;
    delete [] b;
    delete [] Free;
    delete [] Fixed;
    delete [] AphiFixed;
    delete [] buffer;

    MPI_Finalize();

    return 0;
}

void	readData(char* filename, double**& Points, int**& Faces, int**& Elements, Boundary*& Boundaries, int& myN_p, int& myN_f, int& myN_e, int& myN_b, bool*& yourPoints, int myID)
{
    fstream		file;
    string      temp;
    char		myFileName[64]; //Added
    int			myMaxN_sp	= 0;//Added
    int			myMaxN_sb	= 0;//Added
    int			maxN_sp		= 0;//Added
    int			yourID		= 0;//Added

    if(myID==0)
	{
		cout << "Reading " << filename << "'s... " << flush;
	}

    sprintf(myFileName, "%s%d", filename, myID);
    file.open(myFileName);

    file >> temp >> myN_p;
	file >> temp >> myN_f;
	file >> temp >> myN_e;
	file >> temp >> myN_b;

	Points			= new double*	[myN_p];
	Faces			= new int*		[myN_f];
	Elements		= new int*		[myN_e];
	Boundaries		= new Boundary	[myN_b];
    Points[0]       = new double	[myN_p*3];
    Faces[0]        = new int		[myN_f*3];
    Elements[0]     = new int		[myN_e*4];
    yourPoints		= new bool		[myN_p];
    for(int p=1, pp=3; p<myN_p; p++, pp+=3)
    {
        Points[p] = &Points[0][pp];
    }
    for(int f=1, ff=3; f<myN_f; f++, ff+=3)
    {
        Faces[f] = &Faces[0][ff];
    }
    for(int e=1, ee=4; e<myN_e; e++, ee+=4)
    {
        Elements[e] = &Elements[0][ee];
    }
    memset(yourPoints, false, myN_p*sizeof(bool)); //Added

    file >> temp;
    for(int p=0; p<myN_p; p++)
    {
        file >> Points[p][0] >> Points[p][1] >> Points[p][2];
    }

    file >> temp;
    for(int f=0; f<myN_f; f++)
    {
        file >> Faces[f][0] >> Faces[f][1] >> Faces[f][2];
    }

    file >> temp;
    for(int e=0; e<myN_e; e++)
    {
        file >> Elements[e][0] >> Elements[e][1] >> Elements[e][2] >> Elements[e][3]; //Change it! added: Elements[e][3]
    }

    file >> temp;
    for(int b=0; b<myN_b; b++)
	{
		file >> Boundaries[b].name_ >> Boundaries[b].type_ >> Boundaries[b].N_;
		Boundaries[b].indices_  = new int [Boundaries[b].N_];
		for(int n=0; n<Boundaries[b].N_; n++)
		{
			file >> Boundaries[b].indices_[n];
		}
		file >> Boundaries[b].value_;
		if(Boundaries[b].type_=="interprocess")
		{
			myMaxN_sb++;
			myMaxN_sp	= max(myMaxN_sp, Boundaries[b].N_);
			yourID		= static_cast<int> (Boundaries[b].value_);
			if(yourID>myID)
			{
				for(int p=0; p<Boundaries[b].N_; p++)
				{
					yourPoints[Boundaries[b].indices_[p]]	= true;
				}
			}
		}
	}

    MPI_Allreduce(&myMaxN_sp, &maxN_sp, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	buffer		= new double [maxN_sp];
	bufferSize	= (maxN_sp*sizeof(double)+MPI_BSEND_OVERHEAD)*myMaxN_sb;
	MPI_Buffer_attach(new char[bufferSize] , bufferSize);

	file.close();

	if(myID==0)
	{
		cout << "Done.\n" << flush;
	}
	return;
}

void	writeData(double* phi, double**& Points, int**& Elements, int& myN_p, int& myN_e, int l, int myID)
{	
	fstream         file;
    char            fileName[128];
    
    sprintf(fileName, "VTKOutput/mpi/Assignment2_MPI_%02d_%04d.vtk", myID, l);

    file.open(fileName, ios::out);
	
    file << "# vtk DataFile Version 2.0"<< endl;
    file << "untitled, Created by me"	<< endl;
    file << "ASCII"						<< endl;
    file << "DATASET UNSTRUCTURED_GRID"	<< endl;
    
    file << "POINTS " << myN_p << " double" << endl;
    for(int p=0; p<myN_p; p++)
    {
        file << setw(6) << setprecision(5) << fixed << Points[p][0] << "\t" << Points[p][1] << "\t" << Points[p][2] << endl;
    }
    
    file << "CELLS " << myN_e << " " << 5*myN_e << endl;
    for(int e=0; e<myN_e; e++)
    {
        file << "4\t" << Elements[e][0] << "\t" << Elements[e][1] << "\t" << Elements[e][2] <<  "\t" << Elements[e][3] << endl;
    }
    
    file << "CELL_TYPES " << myN_e << endl;
    for(int e=0; e<myN_e; e++)
    {
        file << "10" << endl;
    }
	
    file << "POINT_DATA " << myN_p << endl;
    file << "SCALARS phi double 1" << endl;
    file << "LOOKUP_TABLE \"default\"" << endl;
    for(int p=0; p<myN_p; p++)
    {
        file << setprecision(5) << phi[p] << endl;
    }
	
    file.close();
	
	return;
}

void	exchangeData(double* v, Boundary* Boundaries, int myN_b)
{
	int			yourID		= 0;
	int			tag			= 0;
	MPI_Status	status;

	for(int b=0; b<myN_b; b++)
	{
		if(Boundaries[b].type_=="interprocess")
		{
			for(int p=0; p<Boundaries[b].N_; p++)
			{
				buffer[p]	= v[Boundaries[b].indices_[p]];
			}
			yourID		= static_cast<int> (Boundaries[b].value_);
			MPI_Bsend(buffer, Boundaries[b].N_,	MPI_DOUBLE,	yourID,	tag, MPI_COMM_WORLD);
		}
	}
	for(int b=0; b<myN_b; b++)
	{
		if(Boundaries[b].type_=="interprocess")
		{
			yourID	= static_cast<int> (Boundaries[b].value_);
			MPI_Recv(buffer, Boundaries[b].N_,	MPI_DOUBLE,	yourID,	tag, MPI_COMM_WORLD, &status);
			for(int p=0; p<Boundaries[b].N_; p++)
			{
				v[Boundaries[b].indices_[p]] += buffer[p];
			}
		}
	}

	return;
}


void	assembleSystem(SparseMatrix& M, SparseMatrix& K, double* s, double* phi, bool* Free, bool* Fixed, double** Points, int** Faces, int** Elements, Boundary* Boundaries, int myN_p, int myN_f, int myN_e, int myN_b, int myID)
{
	if(myID==0)
	{
		cout << "Assembling system... " << flush;
	}

    double	x[4];
    double	y[4];
    double	z[4]; //Added
    double	gradEta[3][4]; //Resized
    double	gradEta_p[3]= {0.0, 0.0, 0.0}; //Resized
    double	gradEta_q[3]= {0.0, 0.0, 0.0}; //Resized
    double	M_e[4][4]	= {{2.0, 1.0, 1.0, 1.0}, {1.0, 2.0, 1.0, 1.0}, {1.0, 1.0, 2.0, 1.0}, {1.0, 1.0, 1.0, 2.0}}; //Resized
    double	k_e[3][3]	= {{2.0, 1.0, 1.0}, {1.0, 2.0, 1.0}, {1.0, 1.0, 2.0}}; //Added
    double	s_e[4]		= {1.0, 1.0, 1.0, 1.0}; //Resized
    int		Nodes[4]	= {0, 0, 0, 0}; //Resized
    double* Omega		= new double [myN_e];
    double* Gamma		= new double [myN_f];
    int		m;
    int		n;

    // Assign all the indices to be free initially
    for(int p=0; p<myN_p; p++)
    {
        Free[p] = 1;
        Fixed[p]= 0;
        s[p]    = 0.0;
    }

    // Calculate face areas (Resized)
    for(int f=0; f<myN_f; f++)
    {
        for(int p=0; p<3; p++)
        {
            x[p]	= Points[Faces[f][p]][0];
            y[p]	= Points[Faces[f][p]][1];
            z[p]	= Points[Faces[f][p]][2];
        }
        Gamma[f]	= sqrt(pow(((y[1]-y[0])*(z[2]-z[0]) - (z[1]-z[0])*(y[2]-y[0])),2)
        					+ pow(((z[1]-z[0])*(x[2]-x[0]) - (x[1]-x[0])*(z[2]-z[0])),2)
        					+ pow(((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0])),2))/2;
    }

    // Calculate element volumes (Resized)
    for(int e=0; e<myN_e; e++)
    {
        for(int p=0; p<4; p++)
        {
            x[p]	= Points[Elements[e][p]][0];
            y[p]	= Points[Elements[e][p]][1];
            z[p]	= Points[Elements[e][p]][2];
        }
        Omega[e]    = abs(x[0]*y[1]*z[2] - x[0]*y[2]*z[1] - x[1]*y[0]*z[2]
					 + x[1]*y[2]*z[0] + x[2]*y[0]*z[1] - x[2]*y[1]*z[0]
					 - x[0]*y[1]*z[3] + x[0]*y[3]*z[1] + x[1]*y[0]*z[3]
					 - x[1]*y[3]*z[0] - x[3]*y[0]*z[1] + x[3]*y[1]*z[0]
					 + x[0]*y[2]*z[3] - x[0]*y[3]*z[2] - x[2]*y[0]*z[3]
					 + x[2]*y[3]*z[0] + x[3]*y[0]*z[2] - x[3]*y[2]*z[0]
					 - x[1]*y[2]*z[3] + x[1]*y[3]*z[2] + x[2]*y[1]*z[3]
					 - x[2]*y[3]*z[1] - x[3]*y[1]*z[2] + x[3]*y[2]*z[1]) /6;
    }

    // Assemble M, K, and s
    M.initialize(myN_p, 10);
    K.initialize(myN_p, 10);

    for(int e=0; e<myN_e; e++)
    {
        for(int p=0; p<4; p++)
        {
            Nodes[p]= Elements[e][p];
            x[p]	= Points[Nodes[p]][0];
            y[p]	= Points[Nodes[p]][1];
            z[p]	= Points[Nodes[p]][2]; //Resized
        }

        gradEta[0][0] = ((y[3]-y[1])*(z[2]-z[1])-(y[2]-y[1])*(z[3]-z[1]))/(6*Omega[e]);
		gradEta[0][1] = ((y[2]-y[0])*(z[3]-z[2])-(y[2]-y[3])*(z[0]-z[2]))/(6*Omega[e]);
		gradEta[0][2] = ((y[1]-y[3])*(z[0]-z[3])-(y[0]-y[3])*(z[1]-z[3]))/(6*Omega[e]);
		gradEta[0][3] = ((y[0]-y[2])*(z[1]-z[0])-(y[0]-y[1])*(z[2]-z[0]))/(6*Omega[e]);

		gradEta[1][0] = ((x[2]-x[1])*(z[3]-z[1])-(x[3]-x[1])*(z[2]-z[1]))/(6*Omega[e]);
		gradEta[1][1] = ((x[3]-x[2])*(z[2]-z[0])-(x[0]-x[2])*(z[2]-z[3]))/(6*Omega[e]);
		gradEta[1][2] = ((x[0]-x[3])*(z[1]-z[3])-(x[1]-x[3])*(z[0]-z[3]))/(6*Omega[e]);
		gradEta[1][3] = ((x[1]-x[0])*(z[0]-z[2])-(x[2]-x[0])*(z[0]-z[1]))/(6*Omega[e]);

		gradEta[2][0] = ((x[3]-x[1])*(y[2]-y[1])-(x[2]-x[1])*(y[3]-y[1]))/(6*Omega[e]);
		gradEta[2][1] = ((x[2]-x[0])*(y[3]-y[2])-(x[2]-x[3])*(y[0]-y[2]))/(6*Omega[e]);
		gradEta[2][2] = ((x[1]-x[3])*(y[0]-y[3])-(x[0]-x[3])*(y[1]-y[3]))/(6*Omega[e]);
		gradEta[2][3] = ((x[0]-x[2])*(y[1]-y[0])-(x[0]-x[1])*(y[2]-y[0]))/(6*Omega[e]);

        // Outer loop over each node
        for(int p=0; p<4; p++)
        {
            m		= Nodes[p];
            gradEta_p[0]	= gradEta[0][p];
            gradEta_p[1]	= gradEta[1][p];
            gradEta_p[2]	= gradEta[2][p];

            // Inner loop over each node
            for(int q=0; q<4; q++)
            {
                n			= Nodes[q];
                gradEta_q[0]		= gradEta[0][q];
                gradEta_q[1]		= gradEta[1][q];
                gradEta_q[2]		= gradEta[2][q];

                M(m,n)	   += rho*C*M_e[p][q]*Omega[e]/20;
                K(m,n)	   -= k*(gradEta_p[0]*gradEta_q[0]+gradEta_p[1]*gradEta_q[1]+gradEta_p[2]*gradEta_q[2])*Omega[e];
            }
        }
    }

    // Apply boundary conditions
    for(int b=0; b<myN_b; b++)
    {
        if		(Boundaries[b].type_=="neumann")
        {
            for(int f=0; f<Boundaries[b].N_; f++)
            {
                for(int p=0; p<3; p++)
                {
                    Nodes[p]	= Faces[Boundaries[b].indices_[f]][p];
                    m			= Nodes[p];
                    s[m]	   += s_e[p]*Qcpu*Gamma[Boundaries[b].indices_[f]]/3; //second term contribution to the load vector
                }
            }
        }
        else if	(Boundaries[b].type_=="robin")
        {
            for(int f=0; f<Boundaries[b].N_; f++)
            {
            	for(int p=0; p<3; p++)
				{
					Nodes[p]	= Faces[Boundaries[b].indices_[f]][p];
					m			= Nodes[p];
					s[m]	   += s_e[p]*h*Tair*Gamma[Boundaries[b].indices_[f]]/3; //second term contribution to the load vector

					 for (int q=0; q<3; q++)
					 {
						 Nodes[q]	= Faces[Boundaries[b].indices_[f]][q];
						 n = Nodes[q];
						 K(m,n)    -= h*k_e[p][q]*Gamma[Boundaries[b].indices_[f]]/12; //Contribution to the stiffness matrix
					 }

				}
                Free[m] = true;
            }
        }
    }

    K.finalize();
    M.finalize();

    delete [] Gamma;
    delete [] Omega;

    MPI_Barrier(MPI_COMM_WORLD);

	if(myID==0)
	{
		cout << "Done.\n" << flush;
	}

    return;
}

void	solve(SparseMatrix& A, double* phi, double* b, bool* Free, bool* Fixed, Boundary* Boundaries, bool* yourPoints, int myN_b, int myID)
{
    int		N_row			= A.getNrow();
    double*	r_old			= new double [N_row];
    double*	r				= new double [N_row];
    double*	d				= new double [N_row];
    double*	Ad				= new double [N_row];
    double*	Aphi			= new double [N_row];
    double	alpha			= 0.0;
    double	beta			= 0.0;
    double	r_norm			= 0.0;
    double	tolerance		= 1e-8;
    double	maxIterations	= 1e+3;
	double	r_oldTr_old		= 0.0;
	double	rTr				= 0.0;
	double	dTAd			= 0.0;
	int		k				= 0;
	int		m				= 0;
	//int		n				= 0;

    memset(r_old,		0, N_row*sizeof(double));
    memset(r,			0, N_row*sizeof(double));
    memset(d,			0, N_row*sizeof(double));
    memset(Ad,			0, N_row*sizeof(double));

    // Compute the initial residual
	A.multiply(Aphi, phi, Free, Free);
	exchangeData(Aphi, Boundaries, myN_b); //Added
    for(m=0; m<N_row; m++)
    {
        if(Free[m])
        {
            r_old[m]	= b[m] - Aphi[m];
            d[m]		= r_old[m];
        }
    }
    r_oldTr_old	= computeInnerProduct(r_old, r_old, Free, yourPoints, N_row);
    r_norm		= sqrt(r_oldTr_old);

    // Conjugate Gradient iterative loop
    while(r_norm>tolerance && k<maxIterations)
    {
    	A.multiply(Ad, d, Free, Free);
		exchangeData(Ad, Boundaries, myN_b);
		dTAd	= computeInnerProduct(d, Ad, Free, yourPoints, N_row);
		alpha  	= r_oldTr_old/dTAd;
		for(m=0; m<N_row; m++)
		{
			if(Free[m])
			{
				phi[m] += alpha*d[m];
			}
		}
		for(m=0; m<N_row; m++)
		{
			if(Free[m])
			{
				r[m]	= r_old[m] - alpha*Ad[m];
			}
		}
		rTr	= computeInnerProduct(r, r, Free, yourPoints, N_row);
		beta  	= rTr/r_oldTr_old;
		for(m=0; m<N_row; m++)
		{
			if(Free[m])
			{
				d[m] = r[m] + beta*d[m];
			}
		}
		for(m=0; m<N_row; m++)
		{
			if(Free[m])
			{
				r_old[m] = r[m];
			}
		}
		r_oldTr_old	= rTr;
		r_norm		= sqrt(rTr);
		k++;
	}
	if(myID==0)
	{
		cout << ", k = " << k << ", r_norm = " << r_norm << endl;
	}

    delete [] r_old;
    delete [] r;
    delete [] d;
    delete [] Ad;
    delete [] Aphi;

    return;
}

double	computeInnerProduct(double* v1, double* v2, bool* Free, bool* yourPoints, int N_row)
{
	double		myInnerProduct	= 0.0;
	double		innerProduct	= 0.0;

	for(int m=0; m<N_row; m++)
	{
		if(Free[m] && !yourPoints[m])
		{
			myInnerProduct += v1[m]*v2[m];
		}
	}

	MPI_Allreduce(&myInnerProduct, &innerProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return innerProduct;
}















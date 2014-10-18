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
// Compilation:	g++ Assignment2.cpp -o Assignment2
//
// Execution:	./Assignment2 Box.grid
//
/////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "SparseMatrix.h"
#include "Boundary.h"

using namespace std;

// Global variables
const double	t_min		= 0.00;
const double	t_max		= 100.00;
const double    Delta_t		= 0.1;

const double rho            = 8954.00;
const double C              = 380.00;
const double k              = 386.00;
const double h              = 100.00;
const double Tair           = 300.00;
const double Qcpu           = 40000.00;

const int		N_t			= static_cast<int> ((t_max-t_min)/Delta_t+1);

// Function declarations
void	read(char* filename, double**& Points, int**& Faces, int**& Elements, Boundary*& Boundaries, int& N_p, int& N_f, int& N_e, int& N_b);
void	write(fstream& file, double* phi, int N_p);
void	writeData(double* phi, double**& Points, int**& Elements, int& myN_p, int& myN_e, int l);
void	assemble(SparseMatrix& M, SparseMatrix& K, double* s, double* phi, bool* Free, bool* Fixed, double** Points, int** Faces, int** Elements, Boundary* Boundaries, int N_p, int N_f, int N_e, int N_b);
void	solve(SparseMatrix& A, double* phi, double* b, bool* Free, bool* Fixed);

int     main(int argc, char** argv)
{
    // Simulation parameters
    double**		Points		= NULL;
    int**			Faces		= NULL;
    int**			Elements	= NULL;
    Boundary*		Boundaries	= NULL;
    int				N_p			= 0;
    int				N_f			= 0;
    int				N_e			= 0;
    int				N_b			= 0;
    double			t			= 0;
    fstream         file;

    if(argc<2)
    {
        cerr << "No grid file specified" << endl;
        exit(1);
    }
    else
    {
        read(argv[1], Points, Faces, Elements, Boundaries, N_p, N_f, N_e, N_b);
    }
    // Allocate arrays
    double*			phi			= new double[N_p];
    double*			s			= new double[N_p];
    double*			b			= new double[N_p];
    bool*			Free        = new bool	[N_p];
    bool*			Fixed       = new bool	[N_p];
    double*			AphiFixed	= new double[N_p];
    SparseMatrix	M;
    SparseMatrix	K;
    SparseMatrix	A;

    // Set initial condition
    t			= t_min;
    for(int m=0; m<N_p; m++)
    {
        phi[m]	= 300;
    }

    assemble(M, K, s, phi, Free, Fixed, Points, Faces, Elements, Boundaries, N_p, N_f, N_e, N_b);

    A = M;
    A.subtract(Delta_t, K); // At this point we have A = M-Delta_t*K

    // Compute the column vector to subtract from the right hand side to take account of fixed nodes
    A.multiply(AphiFixed, phi, Free, Fixed);

    //file.open("Assignment2.data", ios::out);
    //write(file, phi, N_p);
    writeData(phi, Points, Elements, N_p, N_e, 0);

    // Time marching loop
    for(int l=0; l<N_t-1; l++)
    {
        t	+= Delta_t;
        cout << "t = " << t;

        // Assemble b
        M.multiply(b, phi);
        for(int m=0; m<N_p; m++)
        {
            b[m]	+= Delta_t*s[m] - AphiFixed[m];
        }

        // Solve the linear system
        solve(A, phi, b, Free, Fixed);

        // Write the solution
    //    write(file, phi, N_p);
        if (l%10==0){
           writeData(phi, Points, Elements, N_p, N_e, (l+1));
        }
    }

    //file.close();

    // Deallocate arrays
    for(int boundary=0; boundary<N_b; boundary++)
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

    return 1;
}

void	read(char* filename, double**& Points, int**& Faces, int**& Elements, Boundary*& Boundaries, int& N_p, int& N_f, int& N_e, int& N_b)
{
    fstream		file;
    string      temp;

    cout << "Reading " << filename << "... " << flush;

    file.open(filename);
    if(!file.is_open())
    {
        cerr << "Error openfing file" << endl;
        exit(1);
    }

    file >> temp >> N_p;
    file >> temp >> N_f;
    file >> temp >> N_e;
    file >> temp >> N_b;

    Points			= new double*	[N_p];
    Faces			= new int*		[N_f];
    Elements		= new int*		[N_e];
    Boundaries		= new Boundary  [N_b];
    Points[0]       = new double    [N_p*3];
    Faces[0]        = new int       [N_f*3];
    Elements[0]     = new int       [N_e*4];

    for(int p=1, pp=3; p<N_p; p++, pp+=3)
    {
        Points[p] = &Points[0][pp];
    }

    for(int f=1, ff=3; f<N_f; f++, ff+=3)
    {
        Faces[f] = &Faces[0][ff];
    }

    for(int e=1, ee=4; e<N_e; e++, ee+=4)
    {
        Elements[e] = &Elements[0][ee];
    }

    file >> temp;
    for(int p=0; p<N_p; p++)
    {
        file >> Points[p][0] >> Points[p][1] >> Points[p][2];
    }

    file >> temp;
    for(int f=0; f<N_f; f++)
    {
        file >> Faces[f][0] >> Faces[f][1] >> Faces[f][2];
    }

    file >> temp;
    for(int e=0; e<N_e; e++)
    {
        file >> Elements[e][0] >> Elements[e][1] >> Elements[e][2] >> Elements[e][3]; //Change it! added: Elements[e][3]
    }

    file >> temp;
    for(int b=0; b<N_b; b++)
    {
        file >> Boundaries[b].name_ >> Boundaries[b].type_ >> Boundaries[b].N_;
        Boundaries[b].indices_  = new int [Boundaries[b].N_];
        for(int n=0; n<Boundaries[b].N_; n++)
        {
            file >> Boundaries[b].indices_[n];
        }
        file >> Boundaries[b].value_;
    }

    file.close();

    cout << "Done.\n" << flush;

    return;
}

void	write(fstream& file, double* phi, int N_p)
{
    for(int m=0; m<N_p; m++)
    {
        file << phi[m] << "\t";
    }
    file << "\n";
    return;
}

void	writeData(double* phi, double**& Points, int**& Elements, int& myN_p, int& myN_e, int l)
{	
	fstream         file;
    char            fileName[64];
    
    sprintf(fileName, "VTKOutput/sequential/Assignment2_%04d.vtk", l);
	
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

void	assemble(SparseMatrix& M, SparseMatrix& K, double* s, double* phi, bool* Free, bool* Fixed, double** Points, int** Faces, int** Elements, Boundary* Boundaries, int N_p, int N_f, int N_e, int N_b)
{
    cout << "Assembling system... " << flush;

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
    double* Omega		= new double [N_e];
    double* Gamma		= new double [N_f];
    int		m;
    int		n;

    // Assign all the indices to be free initially
    for(int p=0; p<N_p; p++)
    {
        Free[p] = 1;
        Fixed[p]= 0;
        s[p]    = 0.0;
    }

    // Calculate face areas (Resized)
    for(int f=0; f<N_f; f++)
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
    for(int e=0; e<N_e; e++)
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
    M.initialize(N_p, 10);
    K.initialize(N_p, 10);

    for(int e=0; e<N_e; e++)
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
            //s[m]		   += s_e[p]*psi*Omega[e]/3;
        }
    }

    // Apply boundary conditions
    for(int b=0; b<N_b; b++)
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
                //Fixed[m]= false;
            }
        }
    }

    K.finalize();
    M.finalize();

    delete [] Gamma;
    delete [] Omega;

    cout << "Done.\n" << flush;

    return;
}

void	solve(SparseMatrix& A, double* phi, double* b, bool* Free, bool* Fixed)
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
    double	N_k	= 1e+3;
    double	r_oldTr_old		= 0.0;
    double	rTr				= 0.0;
    double	dTAd			= 0.0;
    int		k				= 0;
    int		m				= 0;
    int		n				= 0;

    memset(r_old,		0, N_row*sizeof(double));
    memset(r,			0, N_row*sizeof(double));
    memset(d,			0, N_row*sizeof(double));
    memset(Ad,			0, N_row*sizeof(double));

    // Compute the initial residual
    A.multiply(Aphi, phi, Free, Free);
    for(m=0; m<N_row; m++)
    {
        if(Free[m])
        {
            r_old[m]	= b[m] - Aphi[m];
            d[m]		= r_old[m];
            r_oldTr_old+= r_old[m]*r_old[m];
        }
    }
    r_norm = sqrt(r_oldTr_old);

    // Conjugate Gradient iterative loop
    while(r_norm>tolerance && k<N_k)
    {
        dTAd	= 0.0;
        A.multiply(Ad, d, Free, Free);
        for(m=0; m<N_row; m++)
        {
            if(Free[m])
            {
                dTAd   += d[m]*Ad[m];
            }
        }
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
        rTr	= 0.0;
        for(m=0; m<N_row; m++)
        {
            if(Free[m])
            {
                rTr  += r[m]*r[m];
            }
        }
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

    cout << ", k = " << k << ", r_norm = " << r_norm << endl;

    delete [] r_old;
    delete [] r;
    delete [] d;
    delete [] Ad;
    delete [] Aphi;

    return;
}


















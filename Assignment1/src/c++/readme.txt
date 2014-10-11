#Sequential Program: Shallow_Water_Sequential.cpp

	Requierements: Directory Assigment1_Files must be in the same path where Shallow_Water_Sequential.cpp is
	Compile: g++ -o Shallow_Water_Sequential Shallow_Water_Sequential.cpp
	Run Local: ./Shallow_Water_Sequential
	Visualize: Open .csv files from Assigment1_Files directory in Paraview 

#openMP Program: Shallow_Water_OMP.cpp

	Requierements: Directory Assigment1_Files_OMP must be in the same path where Shallow_Water_OMP.cpp is
	Compile: g++ -fopenmp -o Shallow_Water_OMP Shallow_Water_OMP.cpp
	Run Local: ./Shallow_Water_OMP
	Schedule job: sbatch Shallow_Water_OMP.sbatch.avoca
	Visualize: Open .csv files from Assigment1_Files_OMP directory in Paraview 

#MPI Program: Shallow_Water_MPI.cpp

	Requierements: Directory Assigment1_Files_MPI must be in the same path where Shallow_Water_MPI.cpp is
	Compile: mpicxx -o Shallow_Water_MPI Shallow_Water_MPI.cpp
	Run Local: mpirun -np 4 Shallow_Water_MPI
	Schedule job: sbatch Shallow_Water_MPI.sbatch.avoca
	Visualize:  ./postProcessMPI.sh
				Then Open .csv files from Assigment1_Files directory in Paraview 



/**
 * @author github.com/toshNaik
**/

/*
Experiment No. 7
Aim- Write a MPI Parallel Program for Finite Difference Methods
[Textbook section 13.3]
Displays matrix containing the displacement of vibrating string at each time and space interval

Name: Ashutosh Naik
UID: 2018130030
Batch: C

To compile and run:
mpic++ -O -o expt7 expt7.cpp
mpirun -n 2 expt7 m n PRINT

* TODO:: The condition kc/h <= 1 is not checked. 
* resulting in errors if that condition is not satisfied [Explanation in section 13.3.1]
*/

#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>

using namespace std;

// Useful macros
#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p*index+1)-1)/n)

#define F(x) sin(3.14159*(x)) 
#define G(x) 0.0 
#define a 1.0 
#define c 1.0 
#define T 1.0 
int m, n, PRINT;
MPI_Datatype dt_temp, dt_column, sp_temp, sp_column;

// Function declarations
void split_columnwise(double split[], int id, int p);
void gather_columnwise_and_print(double split[], int id, int p);
void print_matrix(double matrix[], int m, int n);
void calculate(double split[], int id, int p);

int main(int argc, char* argv[])
{
	int id, p;
	double start_time;

	if(argc < 4)
		exit(1);
	m = atoi(argv[1]);
	n = atoi(argv[2]);
	PRINT = atoi(argv[3]);

	// MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Add Barrier and start timer	
	MPI_Barrier(MPI_COMM_WORLD);
	if(!id) start_time = MPI_Wtime();

	// Created new MPI type required for columnwise scatter and gather
	MPI_Type_vector(m+1, 1, n+1, MPI_DOUBLE, &dt_temp);
	MPI_Type_create_resized(dt_temp, 0, sizeof(double), &dt_column);
	MPI_Type_commit(&dt_column);

	// Calculating single time and space segments 
	double h = a/n;
	double k = T/m;
	double L = (k*c/h)*(k*c/h);

	// Each process' share of the matrix 
	double *split = (double*)malloc((m+1)*BLOCK_SIZE(id,p,n+1)*sizeof(double*));

	// Calculate the matrix and split it columnwise
	split_columnwise(split, id, p);
	// Calculate the values of the matrix
	calculate(split, id, p);
	// Gather the values from all process' and print the result
	gather_columnwise_and_print(split, id, p);

	if(!id)
	{
		cout << "Time for " << p << " process: ";
		cout << (MPI_Wtime() - start_time)*1000 << " ms" << endl;
	}
	// Return from the code
	MPI_Finalize();
	return 0;
}

void split_columnwise(double split[], int id , int p)
{
	int sum = 0;
	int *sendcounts = (int*)malloc(sizeof(int)*p);
	int *displs = (int*)malloc(sizeof(int)*p);
	
	// Created new MPI type required for columwise scatter and gather
	MPI_Type_vector(m+1, 1, BLOCK_SIZE(id,p,n+1), MPI_DOUBLE, &sp_temp);
	MPI_Type_create_resized(sp_temp, 0, sizeof(double), &sp_column);
	MPI_Type_commit(&sp_column);

	// This loop calculates the displacement and count of the various blocks to be sent to other processes
	for(int i = 0; i < p; i++) 
	{
		sendcounts[i] = BLOCK_SIZE(i, p, n+1);
		displs[i] = sum;
		sum += sendcounts[i];
	}
	if(!id)		// If process id is 0 then calculate the array u and scatter it to other process
	{
		double h = a/n;
		double k = T/m;
		double L = (k*c/h)*(k*c/h);

		double *u = (double*)calloc((m+1)*(n+1), sizeof(double*));

		for(int i=1; i<n; i++)
			u[0*(n+1)+i] = F(i*h);
		
		for(int i=1; i<n; i++)
			u[1*(n+1)+i] = (L/2) * (u[0*(n+1)+(i + 1)] + u[0*(n+1)+(i - 1)]) + (1 - L) * (u[0*(n+1)+ i]) + (k * G(i*h));
	
		MPI_Scatterv(u, sendcounts, displs, dt_column, split, sendcounts[0], sp_column, 0, MPI_COMM_WORLD);
		free(u);
	}
	else 			// Any other process 
		MPI_Scatterv(NULL, NULL, NULL, NULL, split, sendcounts[id], sp_column, 0, MPI_COMM_WORLD);
}

void print_matrix(double matrix[], int m, int n)
{
	// Simple matrix printing function
	cout << "      ";
	for(int i = 0; i < n; i++)
	{
		cout << "|" << i <<"|   ";
	}
	cout << endl;
	for(int i = 0; i < m; i++) {
		cout << "|" << i <<"| ";
		for(int j = 0; j < n; j++)
		{
			printf("%6.2lf", matrix[i * n + j]);
		}
		cout << endl;
	}
	cout << endl;
}  

void gather_columnwise_and_print(double split[], int id, int p)
{
	int sum = 0;
	int *sendcounts = (int*)malloc(sizeof(int)*p);
	int *displs = (int*)malloc(sizeof(int)*p);
	double *u = (double*)calloc((m+1)*(n+1), sizeof(double*));

	// Created new MPI type required for columwise scatter and gather
	MPI_Type_vector(m+1, 1, BLOCK_SIZE(id,p,n+1), MPI_DOUBLE, &sp_temp);
	MPI_Type_create_resized(sp_temp, 0, sizeof(double), &sp_column);
	MPI_Type_commit(&sp_column);

	// This loop calculates the displacement and count of the various blocks to be sent to other processes
	for(int i = 0; i < p; i++) 
	{
		sendcounts[i] = BLOCK_SIZE(i, p, n+1);
		displs[i] = sum;
		sum += sendcounts[i];
	}
	if(!id) 	// If process 0 then gather from other processes into self
	{
		MPI_Gatherv(split, sendcounts[id], sp_column, u, sendcounts, displs, dt_column, 0, MPI_COMM_WORLD);
		if(PRINT)
		{
			cout << "The matrix is: " << endl;
			print_matrix(u, m+1, n+1);
		}
	}
	else		// Other processes
		MPI_Gatherv(split, sendcounts[id], sp_column, u, sendcounts, displs, dt_column, 0, MPI_COMM_WORLD);
}

void calculate(double split[], int id, int p) 
{
	// Calculating single time and space segments 
	double h = a/n;
	double k = T/m;
	double L = (k*c/h)*(k*c/h);
	
	int high = BLOCK_HIGH(id,p,n+1), low = BLOCK_LOW(id,p,n+1);
	
	// Variables to hold ghost values of each row
	double ghost1 = 0, ghost2 = 0;

	for(int j=1; j<m; j++) // Iterating through time segments
	{
		MPI_Request request;
		if(id != 0) // if process 0 then cannot send anything to left
			MPI_Isend(&split[j*BLOCK_SIZE(id,p,n+1) + 0], 1, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, &request);
		if(id != p-1)  // if last process then cannot send anything to right
			MPI_Isend(&split[j*BLOCK_SIZE(id,p,n+1) + high-low], 1, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD, &request);

		if(id != 0) // If not process 0 receive ghost value from left (Process 0 has no left neighbour)
			MPI_Recv(&ghost1, 1, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
		if(id != p-1) // If not last process receive ghost value from right (Process p-1 has no right neighbour)
			MPI_Recv(&ghost2, 1, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

		for(int i=0; i<BLOCK_SIZE(id,p,n+1); i++) // Iterating through space segments
		{
			
			if(id == 0 && i == 0) continue;										// If process 0 and first column then values are zero so continue
			if(id == p-1 && i == BLOCK_SIZE(id,p,n+1)-1) continue;		// If process p-1 and last column then values are zero so continue
			
			if(i == 0) 											// If first column then use ghost1 value
				split[(j+1)*BLOCK_SIZE(id,p,n+1) + i] = (2.0*(1.0 - L) * split[j*BLOCK_SIZE(id,p,n+1) + i]) + (L * (ghost1 + split[j*BLOCK_SIZE(id,p,n+1) + i+1])) - (split[(j-1)*BLOCK_SIZE(id,p,n+1) + i]);	
			else if(i == BLOCK_SIZE(id, p, n+1) - 1) // if last column then use ghost 2 value
				split[(j+1)*BLOCK_SIZE(id,p,n+1) + i] = (2.0*(1.0 - L) * split[j*BLOCK_SIZE(id,p,n+1) + i]) + (L * (split[j*BLOCK_SIZE(id,p,n+1) + i-1] + ghost2)) - (split[(j-1)*BLOCK_SIZE(id,p,n+1) + i]);
			else
				split[(j+1)*BLOCK_SIZE(id,p,n+1) + i] = (2.0*(1.0 - L) * split[j*BLOCK_SIZE(id,p,n+1) + i]) + (L * (split[j*BLOCK_SIZE(id,p,n+1) + i-1] + split[j*BLOCK_SIZE(id,p,n+1) + i+1])) - (split[(j-1)*BLOCK_SIZE(id,p,n+1) + i]);
		}
	}
}
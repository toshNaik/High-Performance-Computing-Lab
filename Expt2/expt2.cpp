/**
 * @author github.com/toshNaik
**/

/*
Experiment No. 2
Aim- Write a MPI Parallel Program to find the shortest path using Floyd's algorithm

Name: Ashutosh Naik
UID: 2018130030
Batch: C

The code reads input from file input.txt
input.txt contains the representation of the graph in adjacency matrix form 
with the first line indicating the number of nodes present in the graph
-1 indicates infinity

To compile and run:
mpic++ -O -o expt2 expt2.cpp
mpirun -n 2 expt2 (PRINT)
*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <mpi.h>

using namespace std;

// Useful macros
#define INF 1000000
#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) ((((p)*((index)+1))-1)/(n))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
int PRINT; 

// Functions
void print_matrix(int matrix[], int n);
void read_row_striped_matrix(int matrix[], int n, int id, int p);
void print_row_striped_matrix(int matrix[], int n, int id, int p);
void compute_shortest_paths(int matrix[], int n, int id, int p);

int main(int argc, char* argv[]) 
{
	double start_time;
	int p, id, n;
	int* matrix;
	PRINT = atoi(argv[1]);

	// MPI Initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	// barrier so all processes start together from this point and then record start time
	MPI_Barrier(MPI_COMM_WORLD);
	if(!id) start_time = MPI_Wtime();

	// If root process: read the size of matrix from file and broadcast it to others
	if(!id)
	{
		FILE * ftr;
		ftr = fopen("input.txt", "r");
		if(ftr == NULL)
			exit(1);
		fscanf(ftr, "%d", &n);
		fclose(ftr);
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// Create the matrix and read from file
	matrix = (int*)malloc(n*BLOCK_SIZE(id, p, n) * sizeof(int));
	read_row_striped_matrix(matrix, n, id, p);
	
	// Floyd's algorithm
	compute_shortest_paths(matrix, n, id, p);
	
	// Gather and print the final matrix
	print_row_striped_matrix(matrix, n, id, p);
	free(matrix);
	
	// If root process: print the execution time to output.txt for plots
	if(!id) {
		ofstream fout("output.txt", ios_base::app);
		fout << "Time for " << p << " processe(s)" << " and " << n << " rows: ";
		fout << (MPI_Wtime() - start_time)*1000 << " ms" << endl;
	}
	MPI_Finalize();
	return 0;
}  

void read_row_striped_matrix(int matrix[], int n, int id, int p)
{ 
	// Reads the matrix from a file
	int value;
	int* temp_mat = NULL;
	int sum = 0;

	int *sendcounts = (int*)malloc(sizeof(int)*p);
   int *displs = (int*)malloc(sizeof(int)*p);

	// This loop calculates the displacement and count of the various blocks to be sent to other processes
	for(int i = 0; i < p; i++) 
	{
		sendcounts[i] = n*BLOCK_SIZE(i, p, n);
		displs[i] = sum;
		sum += sendcounts[i];
	}

	// If process 0 then read from file to temp_mat and scatter blocks to other processes 
	if (!id) {
		temp_mat = (int*)malloc(n*n*sizeof(int));
		int num;
		FILE *ftr;
		ftr = fopen("input.txt", "r");
		fscanf(ftr, "%d", &num);
		for(int i = 0; i < n; i++)	{
			for(int j = 0; j < n; j++) {
				if(i == j)
				{
					fscanf(ftr, "%d", &num);
					temp_mat[i*n+j] = 0;
				}
				else
				{
					fscanf(ftr, "%d", &num);
					if(num == -1)
						temp_mat[i*n+j] = INF;
					else
						temp_mat[i*n+j] = num;
				}
			}
		}
		fclose(ftr);
		MPI_Scatterv(temp_mat, sendcounts, displs, MPI_INT, matrix, sendcounts[id], MPI_INT, 0, MPI_COMM_WORLD);

		// Print the input matrix
		if(PRINT)
		{
			cout << "Input Matrix:" << endl;
			print_matrix(temp_mat, n);
		}
		free(temp_mat);
	} 
	// If not process 0 then receive from process 0
	else 
	{
		MPI_Scatterv(NULL, NULL, NULL, NULL, matrix, sendcounts[id], MPI_INT, 0, MPI_COMM_WORLD);
	}
}

void print_matrix(int matrix[], int n)
{
	// Simple matrix printing function
	if(PRINT) {
		cout << "     ";
		for(int i = 0; i < n; i++)
		{
			cout << "|" << i <<"|  ";
		}
		cout << endl;
		for(int i = 0; i < n; i++) {
			cout << "|" << i <<"| ";
			for(int j = 0; j < n; j++)
			{
				if(matrix[i * n + j] == INF)
				{
					cout << "  inf";
				}
				else
				{
					printf("%5d", matrix[i * n + j]);
				}
				
			}
			cout << endl;
		}
		cout << endl;
	}
}  

void print_row_striped_matrix(int matrix[], int n, int id, int p) 
{
	// Gathering all the local matrices into a common matrix and printing the solution
	int* temp_mat = NULL;

	int sum = 0;

	int *sendcounts = (int*)malloc(sizeof(int)*p);
   int *displs = (int*)malloc(sizeof(int)*p);

	// This loop calculates the displacement and count of the various blocks to be received from other processes
	for (int i = 0; i < p; i++) 
	{
		sendcounts[i] = n*BLOCK_SIZE(i, p, n);
		displs[i] = sum;
		sum += sendcounts[i];
	}

	// If process 0 receive from all processes and print the result
	if (!id) {
		temp_mat = (int*)malloc(n*n*sizeof(int));
		MPI_Gatherv(matrix, sendcounts[0], MPI_INT, temp_mat, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		
		if(PRINT)
		{ 
			cout << "The solution is:" << endl;
			print_matrix(temp_mat, n);
		}
		free(temp_mat);
	} 
	// If not process 0 then send to process 0
	else 
	{
		MPI_Gatherv(matrix, sendcounts[id], MPI_INT, temp_mat, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	}
}

void compute_shortest_paths(int matrix[], int n, int id, int p) 
{
	int* tmp = (int*)malloc(n*sizeof(int));
	for(int k = 0; k < n; k++)
	{
		// Get the process to which row k belongs
		int root = BLOCK_OWNER(k, p, n);
		if(id == root)
		{
			// Offset gets to the correct row within the local matrix
			int offset = k - BLOCK_LOW(id, p, n);
			// Get all the elements of that row into matrix tmp
			for (int j = 0; j < n; j++)
				tmp[j] = matrix[offset * n + j];
		}
		// Broadcast that row to all processes
		MPI_Bcast(tmp, n, MPI_INT, root, MPI_COMM_WORLD);

		// Finding the shortest path
		for(int i = 0; i < BLOCK_SIZE(id, p, n); i++)
			for(int j = 0; j < n; j++) 
				matrix[i*n+j] = MIN(matrix[i*n+j], matrix[i*n+k] + tmp[j]);
	}
	free(tmp);
}   
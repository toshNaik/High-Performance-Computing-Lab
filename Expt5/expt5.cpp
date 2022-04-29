/**
 * @author github.com/toshNaik
**/

/*
Experiment No. 5
Aim- Write a MPI Parallel Program for Matrix Multiplication
[Textbook section 11.3]
Rowwise block-striped parallel algorithm
Works only for nxn matrices 

Name: Ashutosh Naik
UID: 2018130030
Batch: C

The code reads matrices from files matrix1.txt and matrix2.txt
the firt line of each file indicates the number of rows and columns 
followed by the matrix

To compile and run:
mpic++ -O -o expt5 expt5.cpp
mpirun -n 2 expt5 (PRINT)
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
int PRINT;

// Function declarations
void read_row_striped_matrix(int matrix[], int rows, int cols, int id, int p, char filename[]);
void multiply_matrices(int matrix1[], int matrix2[], int result[], int n, int id, int p);
void print_matrix(int matrix[], int n);
void print_row_striped_matrix(int matrix[], int n, int id, int p); 

int main(int argc, char* argv[]) {
	int id, p, n;
	double start_time;
	PRINT = atoi(argv[1]);

	// MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
		
	// Add Barrier and start timer	
	MPI_Barrier(MPI_COMM_WORLD);
	if(!id)	start_time = MPI_Wtime();
	
	// If root process: read the size of matrix
	if(!id)
	{
		FILE* ftr;
		ftr = fopen("matrix1.txt", "r");
		fscanf(ftr, "%d", &n);
		fclose(ftr);
	}
	// Broadcast n to other processes.
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate space for local matrices and result.
	int* matrix1 = (int*)malloc(BLOCK_SIZE(id, p, n)*n*sizeof(int));
	int* matrix2 = (int*)malloc(BLOCK_SIZE(id, p, n)*n*sizeof(int));
	int* result = (int*)malloc(BLOCK_SIZE(id, p, n)*n*sizeof(int));
	for(int i=0; i<BLOCK_SIZE(id, p, n); i++)
		for(int j=0; j<n; j++)
			result[i*n+j] = 0; 

	// Read the matrices from file
	read_row_striped_matrix(matrix1, n, n, id, p, "matrix1.txt");
	read_row_striped_matrix(matrix2, n, n, id, p, "matrix2.txt");

	// Multiply the matrices
	multiply_matrices(matrix1, matrix2, result, n, id, p);

	// Print the resultant matrix
	print_row_striped_matrix(result, n, id, p);

	// If root process: print the execution time to output.txt
	if(!id) {
		ofstream fout("output.txt", ios_base::app);
		fout << "Time for " << p << " process" << " and " << n << " rows: ";
		fout << (MPI_Wtime() - start_time)*1000 << " ms" << endl;
	}

	MPI_Finalize();
	return 0;
}


void read_row_striped_matrix(int matrix[], int rows, int cols, int id, int p, char filename[])
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
		sendcounts[i] = cols*BLOCK_SIZE(i, p, rows);
		displs[i] = sum;
		sum += sendcounts[i];
	}

	// If process 0 then read from file to temp_mat and scatter blocks to other processes 
	if (!id) {
		temp_mat = (int*)malloc(rows*cols*sizeof(int));
		int num;
		FILE *ftr;
		ftr = fopen(filename, "r");
		fscanf(ftr, "%d", &num);
		fscanf(ftr, "%d", &num);
		for(int i = 0; i < rows; i++)	{
			for(int j = 0; j < cols; j++) {
				fscanf(ftr, "%d", &num);
				temp_mat[i*cols+j] = num;
			}
		}
		fclose(ftr);
		MPI_Scatterv(temp_mat, sendcounts, displs, MPI_INT, matrix, sendcounts[id], MPI_INT, 0, MPI_COMM_WORLD);
		free(temp_mat);
	} 
	// If not process 0 then receive from process 0
	else 
	{
		MPI_Scatterv(NULL, NULL, NULL, NULL, matrix, sendcounts[id], MPI_INT, 0, MPI_COMM_WORLD);
	}
}


void multiply_matrices(int matrix1[], int matrix2[], int result[], int n, int id, int p) {
	int block_size_1 = BLOCK_SIZE(id, p, n);
	int block_size_2 = BLOCK_SIZE(id, p, n);
	int block_low_2 = BLOCK_LOW(id, p, n);
	int next_block_size_2, next_block_low_2;

	for(int i=0; i<p; i++) {
		MPI_Request request, request2, request3;

		// Sending the block size and block low to previous process and receiving from next
		MPI_Isend(&block_size_2, 1, MPI_INT, (id-1+p)%p, 35, MPI_COMM_WORLD, &request);
		MPI_Recv(&next_block_size_2, 1, MPI_INT, (id+1)%p, 35, MPI_COMM_WORLD, NULL);
		MPI_Isend(&block_low_2, 1, MPI_INT, (id-1+p)%p, 36, MPI_COMM_WORLD, &request2);
		MPI_Recv(&next_block_low_2, 1, MPI_INT, (id+1)%p, 36, MPI_COMM_WORLD, NULL);

		// Initializing temp buffer to hold matrix2 of next process
		int* temp_2 = (int*)malloc(next_block_size_2*n*sizeof(int));

		// Sending matrix2 to previous process and receiving from next
		MPI_Isend(matrix2, n*block_size_2, MPI_INT, (id-1+p)%p, 37, MPI_COMM_WORLD, &request3);
		MPI_Recv(temp_2, n*next_block_size_2, MPI_INT, (id+1)%p, 37, MPI_COMM_WORLD, NULL);

		// Wait for non blocking send to finish
		MPI_Wait(&request, NULL);
		MPI_Wait(&request2, NULL);
		MPI_Wait(&request3, NULL);

		// Current block length, block low and matrix are received values
		block_size_2 = next_block_size_2;
		block_low_2 = next_block_low_2;
		free(matrix2);
		matrix2 = temp_2;
		
		// The result matrix 
		for(int i=0; i<block_size_1; i++)
			for(int j=block_low_2; j<block_low_2+block_size_2; j++)
				for(int k=0; k<n; k++)
					result[i*n+k] += matrix1[i*n+j] * matrix2[(j-block_low_2)*n+k];
	}
}

void print_matrix(int matrix[], int n)
{
	// Simple matrix printing function
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
			printf("%5d", matrix[i * n + j]);
		}
		cout << endl;
	}
	cout << endl;
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

		if(PRINT) {
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

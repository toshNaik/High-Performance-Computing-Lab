/**
 * @author github.com/toshNaik
**/

/*
Experiment No. 6
Aim- Write a MPI Parallel Program for solving Linear System of Equations
[Textbook section 12.4.3]
Row-oriented algorithm

Name: Ashutosh Naik
UID: 2018130030
Batch: C

The code reads coeffecients from file matrix.txt
the firt line of each file indicates the number of equations
followed by equations of the form: Ax1 + Bx2 + Cx3 + .... = Z

To compile and run:
mpic++ -O -o expt6 expt6.cpp
mpirun -n 2 expt6 (PRINT)
*/

#include<mpi.h>
#include<bits/stdc++.h>

using namespace std;

// Useful macros
#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) ((((p)*((index)+1))-1)/n)
int PRINT;

// Required functions
void read_augmented_matrix(double aug_mat[], int n, int id, int p);
void gaussian_elimination(double aug_mat[], int loc[], int n, int id, int p);
void back_substitution(double aug_mat[], double x[], int loc[], int n, int id, int p);
void gather_result(double result[], int id, int p, int n);

// Struct required for All_reduce operation
typedef struct ind_val
{
	double value;
	int index;
} Pair;

int main(int argc, char* argv[])
{
	int id, p, n;
	double start_time, end_time;
	PRINT = atoi(argv[1]);

	// MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
		
	// Add Barrier and start timer	
	MPI_Barrier(MPI_COMM_WORLD);
	if(!id)	start_time = MPI_Wtime();

	// open file and get number of equations
	FILE* input = fopen("matrix.txt", "r");
	if(!id)
		fscanf(input, "%d", &n);
	fclose(input);

	// Broadcast number of equations to all processes
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	double* aug_mat = (double*)malloc(BLOCK_SIZE(id, p, n)*(n+1)*sizeof(double));
	double* x = (double*)malloc(BLOCK_SIZE(id, p, n)*sizeof(double));
	
	// Read the augmented matrix from file
	read_augmented_matrix(aug_mat, n, id, p);

	// The permutation of the columns of the current row
	int loc[n];
	for(int i=0; i<n; i++)
		loc[i] = i;

	// Gaussian elimination
	gaussian_elimination(aug_mat, loc, n, id, p);
	// Back substitution and gather results
	back_substitution(aug_mat, x, loc, n, id, p);
	gather_result(x, id, p, n);

	end_time = MPI_Wtime() - start_time;

	if(!id) {
		ofstream fout("expt6_output.txt", ios_base::app);
		fout << "Time for " << p << " process" << " and " << n << " rows: ";
		fout << (MPI_Wtime() - start_time)*1000 << " ms" << endl;
	}
	MPI_Finalize();
	return 0;
}

void back_substitution(double aug_mat[], double x[], int loc[], int n, int id, int p)
{
	for(int i=n-1; i>=0; i--)
	{
		int root = BLOCK_OWNER(loc[i], p, n);
		double temp;
		if(id==root)
			temp = aug_mat[(loc[i]-BLOCK_LOW(id,p,n))*(n+1)+n]/aug_mat[(loc[i]-BLOCK_LOW(id,p,n))*(n+1)+i];

		MPI_Bcast(&temp, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

		int row_root = BLOCK_OWNER(i, p, n); // get process actually responsible for holding xi
		
		if(id==row_root){
			x[i-BLOCK_LOW(id,p,n)] = temp; // process actually responsible for holding xi stores its value
		}

		for(int j=0; j<i; j++){ // all processes substitute the value of xi in its equations
			int row_root2 = BLOCK_OWNER(loc[j], p, n); 
			if(id==row_root2){//if current process owns the loc[j] row
				aug_mat[(loc[j]-BLOCK_LOW(id,p,n))*(n+1)+n] -= temp*aug_mat[(loc[j]-BLOCK_LOW(id,p,n))*(n+1)+i];
			}
		}
	}
}

void gaussian_elimination(double aug_mat[], int loc[], int n, int id, int p)
{
	int* marked = (int*)calloc(BLOCK_SIZE(id,p,n), sizeof(int)); // to keep track of marked rows
	// Iterate through each column
	for(int i=0; i<n-1; i++)
	{
		Pair local, global;
		local.value = 0;
		local.index = 0;
		// Iterate through each row 
		for(int j=0; j<BLOCK_SIZE(id, p, n); j++)
		{
			// If row was used as pivot in the past then skip
			if(marked[j]) continue;
			if(abs(aug_mat[j*(n+1)+i]) > local.value)
			{
				local.value = aug_mat[j*(n+1)+i];
				local.index = BLOCK_LOW(id, p, n) + j; // added with BLOCK_LOW to get global index of that row
			}
		}
		// After this op all processes have largest value in i-th column and its global index
		MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
		// find the block owner to mark that row
		int root_marked = BLOCK_OWNER(global.index, p, n); 
		if(root_marked == id)
			marked[global.index-BLOCK_LOW(id,p,n)] = 1;
		
		// adjust loc array to hold the correct order of rows
		int temp = loc[i];
		loc[i] = loc[global.index];
		loc[global.index] = temp;

		double* buffer = (double*)malloc(sizeof(double)*(n+1));
		// Broadcast the pivot row
		int root = BLOCK_OWNER(loc[i], p, n);
		if(id == root)
			for(int j=0; j<n+1; j++)
				buffer[j] = aug_mat[(loc[i]-BLOCK_LOW(id,p,n))*(n+1)+j];
		MPI_Bcast(buffer, n+1, MPI_DOUBLE, root, MPI_COMM_WORLD);

		// drive to 0 the ith column of all unmarked rows in each process
		for(int j=i+1; j<n; j++)
		{
			if(BLOCK_LOW(id,p,n)<=loc[j] && BLOCK_HIGH(id,p,n)>=loc[j])
			{
				double t = aug_mat[(loc[j]-BLOCK_LOW(id,p,n))*(n+1)+i]/buffer[i];
				for(int k=i; k<n+1; k++)
				{
					aug_mat[(loc[j]-BLOCK_LOW(id,p,n))*(n+1)+k] -= buffer[k]*t;
				}
			}
		}
	}
}

void read_augmented_matrix(double aug_mat[], int n, int id, int p)
{ 
	// Reads the matrix from a file
	int value;
	int sum = 0;

	int *sendcounts = (int*)malloc(sizeof(int)*p);
   int *displs = (int*)malloc(sizeof(int)*p);

	// This loop calculates the displacement and count of the various blocks to be sent to other processes
	for(int i = 0; i < p; i++) 
	{
		sendcounts[i] = (n+1)*BLOCK_SIZE(i, p, n);
		displs[i] = sum;
		sum += sendcounts[i];
	}

	// If process 0 then read from file to temp_mat and scatter blocks to other processes 
	if (!id) {
		double* temp_mat = (double*)malloc(n*(n+1)*sizeof(double));
		double num;
		FILE *ftr;
		ftr = fopen("matrix.txt", "r");
		fscanf(ftr, "%lf", &num);
		for(int i = 0; i < n; i++)	{
			// n+1 in inner loop because augmented matrix
			for(int j = 0; j < n+1; j++) {
				fscanf(ftr, "%lf", &num);
				temp_mat[i*(n+1)+j] = num;
			}
		}
		fclose(ftr);

		// Scatter blocks of temp_mat to other processes
		MPI_Scatterv(temp_mat, sendcounts, displs, MPI_DOUBLE, aug_mat, sendcounts[id], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		free(temp_mat);
	} 
	// If not process 0 then receive from process 0
	else 
		MPI_Scatterv(NULL, NULL, NULL, NULL, aug_mat, sendcounts[id], MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void gather_result(double result[], int id, int p, int n)
{
	int sum = 0;
	double* temp_mat = (double*)malloc(n*sizeof(double));
	int* sendcounts = (int*)malloc(p*sizeof(int));
	int* displs = (int*)malloc(p*sizeof(int));
	
	// This loop calculates the displacement and count of the various blocks to be received
	for(int i=0; i<p; i++)
	{
		sendcounts[i] = BLOCK_SIZE(i, p, n);
		displs[i] = sum;
		sum += sendcounts[i];
	}
	MPI_Gatherv(result, sendcounts[id], MPI_DOUBLE, temp_mat, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(PRINT) {
		// Printing the resultant vector
		if(!id) {
			cout << "Result: " << endl;
			for(int i=0; i<n; i++)
				cout << "x" << (i+1) << "\t=\t" << temp_mat[i] << endl;
		}
	}
}
/**
 * @author github.com/toshNaik
**/

/*
Experiment No. 1
Aim- Write a MPI Parallel Program to find all primes numbers less than N where N is large number [Use Sieve of Eratosthenes]

Name: Ashutosh Naik
UID: 2018130030
Batch: C

To compile and run:
mpic++ -O -o expt1 expt1.cpp
mpirun -n 2 expt1 (range)
*/

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p*index+1)-1)/n)

int main(int argc, char *argv[])
{
	int i, id, p, n;
	int high_value, low_value, size;
	int first, prime, index;
	char* marked;
	int global_count, count;
	double start_time, elapsed_time;
	
	// MPI Initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Setup barrier and start timer
	MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime();

	// If no range provided then exit else convert range to int
	if(argc != 2) {
		if(!id) printf("Give a range\n");
		MPI_Finalize();
		exit(1);
	}
	n = atoi(argv[1]);

	// Condition checks if all unmarked primes would be in process 0 
	// [Textbook section 5.4.5]
	// if process 0 would not contain all the unmarked primes then exit
	if((2 + BLOCK_SIZE(0,p,n-1)) < ((int)sqrt((double)n))) {
		if(!id) printf("Error");
		MPI_Finalize();
		exit(1);
	}
	
	// Find the lower and upper values in the block and its size.
	low_value = 2 + BLOCK_LOW(id, p, n-1);
	high_value = 2 + BLOCK_HIGH(id, p, n-1);
	size = BLOCK_SIZE(id, p, n-1);
	
	// Create list of natural numbers for each process
	marked = (char *) malloc(size);
	if (marked == NULL) {
		printf("Cannot allocate enough memory.\n");
		MPI_Finalize();
		exit(1);
	}
	for(i=0; i<size; i++) marked[i] = 0;
	
	// Set k=2 and initialize index of k in process 0
	prime = 2;
	if(!id) index = 0;
	
	// Repeat while k*k <= n
	do {
		// Finding the first element in the block that is a multiple of k
		if(prime * prime > low_value) {
			first = prime * prime - low_value;
		}
		else {
			if(!(low_value % prime)) first = 0;
			else first = prime - (low_value % prime);
		}
		// Marking all the multiples of k in that block
		for(i = first; i < size; i += prime) marked[i] = 1;

		//  If process 0 then find the next smallest unmarked prime and set that as the next k
		if(!id) {
			while (marked[++index]);
			prime = index + 2;
		}

		// Broadcast the value of k found by process 0 to other processes.
		MPI_Bcast(&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
	} while (prime * prime <= n);
	
	// Find the count of  unmarked numbers in each process and then sum it in the reduction step.
	count = 0;
	for(i = 0; i < size; i++)
		if(!marked[i]) count++;
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  	// Stop the timer
	elapsed_time = MPI_Wtime() - start_time;

	// If process 0 then printf to stdout.
	if(!id) {
		printf("In %f seconds we found %d primes less than or equal to %d.\n",
		elapsed_time, global_count, n);
	}

	MPI_Finalize();
	return 0;
}

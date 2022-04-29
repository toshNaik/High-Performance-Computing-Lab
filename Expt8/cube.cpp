/*
	Trying out MPI_Cart functions to create hypercube topology
*/

#include<mpi.h>
#include<bits/stdc++.h>

using namespace std;

int main(int argc, char* argv[])
{
	int rank, size; //I am process RANK and we are a total of SIZE
	MPI_Init(&argc, &argv); 

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	
	MPI_Comm nthCube;
	int nDim=3;
	int processPerDim [3]= {2,2,2};
	int period [3]= {1,1,1};

	MPI_Cart_create(MPI_COMM_WORLD, nDim, processPerDim, period, true, &nthCube);

	int rankInDim;
	MPI_Comm_rank(nthCube, &rankInDim);

	int rank_source, rank_desta, rank_destb, rank_destc, rank_destd;
	MPI_Cart_shift(nthCube, 0,1,&rank_source, &rank_desta);
	MPI_Cart_shift(nthCube, 1,1,&rank_source, &rank_destb);
	MPI_Cart_shift(nthCube, 2,1,&rank_source, &rank_destc);
	for(int i=0; i<size; i++)
	{
		MPI_Barrier(nthCube);
		if(i == rank)
			cerr << "I am known in the world as " << rankInDim << " my adjacents are -> " << rank_desta << " - " << rank_destb << " - " << rank_destc << endl; // << "-" << rank_destd <<"\n";
	}

	// MPI_Comm new_comm1;
	// int remain_dims[] = {0, 1, 1};
	// MPI_Cart_sub(nthCube, remain_dims, &new_comm1);

	// int rank_in;
	// MPI_Comm_rank(new_comm1, &rank_in);
	// MPI_Cart_shift(new_comm1, 0,1,&rank_source, &rank_desta);
	// for(int i=0; i<size; i++)
	// {
		// MPI_Barrier(new_comm1);
		// if(i == rank)
			// cerr << rankInDim << "  " << "I am known in the world as " << rank_in << " my adjacents are -> " << rank_desta << " - " << rank_destb << " - " << rank_destc << endl; // << "-" << rank_destd <<"\n";
	// }
	MPI_Finalize();
	return 0;
}
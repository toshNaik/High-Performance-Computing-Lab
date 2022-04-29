/**
 * @author github.com/toshNaik
**/

/*
Experiment No. 8
Aim- Write a MPI Parallel Program for Hyperquicksort
[Textbook section 14.4]

Name: Ashutosh Naik
UID: 2018130030
Batch: C

Reads numbers from input.txt and sorts them in descending order
The first line of input.txt is the size of array to be sorted
following line contains space separated array

To compile and run:
mpic++ -O -o expt8 expt8.cpp
mpirun -n 4 expt8
*/

#include<mpi.h>
#include<bits/stdc++.h>

using namespace std;

// Useful macros
#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) ((((p)*((index)+1))-1)/n)

// Functions used 
void read_numbers_from_file(vector<int>& portion, string filename, int id, int p, int n);
int compare(const void* a, const void* b);
void hyperquicksort(vector<int>& portion, int id, int p, int n_dims, MPI_Comm comm);
void print_result(vector<int>& portion, int id, int p, MPI_Comm comm);

// comparator function to reverse merge sort
struct greaters {
   bool operator()(const int& a, const int& b) const
   {
		return a > b;
   }
};

int main(int argc, char* argv[])
{
	int p, id, id_cube, n;
	double start_time, end_time;

	// MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// If number of processes is not power of 2 then abort
	if(ceil(log2(p)) != floor(log2(p)))
	{
		cout << "Number of processes must be power of 2" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	// Add Barrier and start timer	
	MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime();

	// Get size of array to be sorted from input.txt
	if(!id)
	{
		ifstream input_file("input.txt");
		input_file >> n;
		input_file.close();
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast it to all processes

	vector<int> portion(BLOCK_SIZE(id,p,n)); // Create portion array to hold portion of numbers
	read_numbers_from_file(portion, "input.txt", id, p, n);

	// Creating communicator for hypercube topology
	MPI_Comm Hypercube_comm;
	int n_dim = (int)log2(p);
	int process_per_dim[n_dim];
	int periods[n_dim];
	for(int i=0; i<n_dim; i++)
	{
		process_per_dim[i] = 2;
		periods[i] = 1;
	}
	MPI_Cart_create(MPI_COMM_WORLD, n_dim, process_per_dim, periods, true, &Hypercube_comm);

	// Get new id in hypercube communicator
	MPI_Comm_rank(Hypercube_comm, &id_cube);
	
	// Print each process' share of array before sorting
	if(!id)
		cout << "\n--- Before sorting --- " << endl;
	print_result(portion, id_cube, p, Hypercube_comm);

	// Sort the portion of data in each process and call hyperquicksort 
	qsort(portion.data(), portion.size(), sizeof(int), compare);
	hyperquicksort(portion, id_cube, p, n_dim, Hypercube_comm);

	// Print each process' share of array after sorting
	if(!id)
		cout << "\n--- After sorting --- " << endl;
	print_result(portion, id_cube, p, Hypercube_comm);

	end_time = MPI_Wtime() - start_time;

	if(!id) {
		ofstream fout("expt8_output.txt", ios_base::app);
		fout << "Time for " << p << " process(es)" << " and " << n << " rows: ";
		fout << (MPI_Wtime() - start_time)*1000 << " ms" << endl;
	}

	MPI_Finalize();
	return 0;
}

void hyperquicksort(vector<int>& portion, int id, int p, int n_dims, MPI_Comm comm)
{
	// Finding the median of root and broadcasting to other processes
	int median;
	if(!id)
		(portion.size()%2 == 0)? median = portion[portion.size()/2 - 1] : median = portion[portion.size()/2];
	MPI_Bcast(&median, 1, MPI_INT, 0, comm);
	
	// Find the position to split the array
	vector<int>::iterator it;
	for(it=portion.begin(); it!=portion.end(); it++)
		if(*it < median)
			break;
	// Split the array into two halves each containing values less than median and greater than median respectively
	vector<int> greater_half(portion.begin(), it);
	vector<int> lesser_half(it, portion.end());

	// Get the ranks of process pairs along a dimension in the hypercube
	int rank_source, rank_dest;
	MPI_Cart_shift(comm, 0, 1, &rank_source, &rank_dest);

	// If rank of destination is higher, send lesser half 
	if(id < rank_dest)
	{
		// Sending size of lesser half to neighbour and receiving size of greater half from neighbour
		int buffer_size;
		int size = lesser_half.size();
		MPI_Sendrecv(&size, 1, MPI_INT, rank_dest, 0, 
						 &buffer_size, 1, MPI_INT, rank_dest, 1, comm, NULL);
		// Sending lesser half to neighbour and receiving greater half from neighbour into buffer
		vector<int> buffer(buffer_size);
		MPI_Sendrecv(&lesser_half[0], size, MPI_INT, rank_dest, 0,
						 &buffer[0], buffer_size, MPI_INT, rank_dest, 1, comm, NULL);

		// Merge sort the two halves 
		portion.resize(greater_half.size() + buffer_size);
		merge(buffer.begin(), buffer.end(), greater_half.begin(), greater_half.end(), portion.begin(), greaters());
	}
	// else rank of destination is lesser, send greater half
	else
	{
		// Sending size of greater half to neighbour and receiving size of lesser half from neighbour
		int buffer_size;
		int size = greater_half.size();
		MPI_Sendrecv(&size, 1, MPI_INT, rank_dest, 1, 
						 &buffer_size, 1, MPI_INT, rank_dest, 0, comm, NULL);
		// Sending greater half to neighbour and receiving lesser half from neighbour into buffer
		vector<int> buffer(buffer_size);
		MPI_Sendrecv(&greater_half[0], size, MPI_INT, rank_dest, 1, 
						 &buffer[0], buffer_size, MPI_INT, rank_dest, 0, comm, NULL);

		// Merge sort the two halves 
		portion.resize(lesser_half.size() + buffer_size);
		merge(buffer.begin(), buffer.end(), lesser_half.begin(), lesser_half.end(), portion.begin(), greaters());
	}

	// Base condition to end recursion
	if(!(n_dims-1))
		return;

	// Partition the hypercube
	int remain_dims[n_dims];
	remain_dims[0] = 0;
	for(int i=1; i<n_dims; i++)
		remain_dims[i] = 1;
	MPI_Comm sub_comm;
	MPI_Cart_sub(comm, remain_dims, &sub_comm);

	// Get the new ranks of the partitioned hypercube
	int sub_id;
	MPI_Comm_rank(sub_comm, &sub_id);

	// Call hyperquicksort on the partitioned hypercubes
	hyperquicksort(portion, sub_id, p/2, n_dims-1, sub_comm);
}

void read_numbers_from_file(vector<int>& portion, string filename, int id, int p, int n)
{
	int* sendcounts = (int*)malloc(p*sizeof(int));
	int* displs = (int*)malloc(p*sizeof(int));
	// This loop calculates the displacement and count of the various blocks to be received
	int sum = 0;
	for(int i=0; i<p; i++)
	{
		sendcounts[i] = BLOCK_SIZE(i, p, n);
		displs[i] = sum;
		sum += sendcounts[i];
	}
	// if root then read from file and scatter 
	if(!id)
	{
		vector<int> array(n);
		// Read from input file the array and store it into array
		ifstream input_file(filename);
		int number;
		input_file >> number;
		for(int i=0; input_file >> number; i++)
			array[i] = number;
		input_file.close();

		MPI_Scatterv(&array[0], sendcounts, displs, MPI_INT, &portion[0], sendcounts[id], MPI_INT, 0, MPI_COMM_WORLD);
	}
	else
		MPI_Scatterv(NULL, NULL, NULL, NULL, &portion[0], sendcounts[id], MPI_INT, 0, MPI_COMM_WORLD);
	free(sendcounts);
	free(displs);
}

void print_result(vector<int>& portion, int id, int p, MPI_Comm comm)
{
	for(int i=0; i<p; i++)
	{
		MPI_Barrier(comm);
		if(i == id)
		{
			cout << "Process " << id << " : ";
			for(auto i: portion)
				cout << i << " ";
			cout << endl;
		}
	}
}

int compare(const void* a, const void* b)
{
	const int* x = (int*) a;
	const int* y = (int*) b;

	return (*y - *x);
}
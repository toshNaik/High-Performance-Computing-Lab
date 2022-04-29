/*
	Sample file to test serialization and deserialization of unordered map 
*/

#include <mpi.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <bits/stdc++.h>

namespace archive = boost::archive;

int main(int argc, char* argv[])
{
	int p, id;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if(!id) {
		std::unordered_map<std::string, int> umap;
		umap.insert({"hello", 1});
		umap.insert({"there", 2});
		umap.insert({"ashutosh", 3});

		std::ostringstream oss;
		archive::binary_oarchive oa(oss);
		oa << umap;
		
		MPI_Send(oss.str().c_str(), oss.str().size(), MPI_BYTE, 1, 32, MPI_COMM_WORLD);
		std::cout << "Map sent by process rank: " << id << std::endl;
	}
	else {
		int size;
		MPI_Status status;

		std::vector<char> receiveBuffer(200);
		MPI_Recv(&receiveBuffer[0], receiveBuffer.size(), MPI_BYTE, 0, 32, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_BYTE, &size);
		receiveBuffer.resize(size);
		std::istringstream iss(std::string(&receiveBuffer[0], receiveBuffer.size()));
		archive::binary_iarchive ia(iss);
		std::unordered_map<std::string, int> umap;
		ia >> umap;
		for(auto i: umap) {
			std::cout << i.first << " " << i.second << std::endl;
		}
	}

	MPI_Finalize();
	exit(0);
}
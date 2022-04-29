/**
 * @author github.com/toshNaik
**/

/*
Experiment No. 4
Aim- Write a MPI Parallel Program for Document Classification

Name: Ashutosh Naik
UID: 2018130030
Batch: C

The code iterates through the directory structure finding .txt files
then prints the profile vector of each txt file (the count of words present in the dictionary)

To compile and run:
mpic++ -O -o expt4 expt4.cpp -lboost_system -lboost_filesystem -lboost_serialization
mpirun --oversubscribe -n 3 expt4 (DIR) (dictionary.txt)
*/

#include <mpi.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/filesystem.hpp>
#include <bits/stdc++.h>

#define DICT_SIZE_MSG 0
#define FILE_NAME_MSG 1
#define VECTOR_MSG 2
#define EMPTY_MSG 3

#define DIR_ARG 1
#define DICT_ARG 2

using namespace std;
namespace archive = boost::archive;
namespace fs = boost::filesystem;

// Function declarations
vector<fs::path> get_names(fs::path const & root, string const & ext);
void make_profile(char* name, unordered_map<string, int>& dict);
void read_dictionary(char* filename, vector<char*>& buffer, int& file_len);

int main(int argc, char* argv[])
{
    int p, id;
    MPI_Comm worker_comm;

    void manager(int, char **, int);
    void worker(int, char **, MPI_Comm);
    
    // MPI Initialization
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

    // If args less than three then exit
    if(argc < 3) {
        cout << "Atleast 3 arguments need to be provided" << endl;
        exit(1);
    }
    // Atleast 2 processes required for worker manager
    if(p < 2) {
        cout << "Atleast 2 processes needed" << endl;
    }
    // Split the comm for global and worker
    else {
        if(!id) {
            MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, id, &worker_comm);
            manager(argc, argv, p);
        }
        else {
            MPI_Comm_split(MPI_COMM_WORLD, 0, id, &worker_comm);
            worker(argc, argv, worker_comm);
        }
    }

    MPI_Finalize();
    return 0;
}	

void manager(int argc, char *argv[], int p) {
    MPI_Request pending;
    MPI_Status status;
    int dict_size, file_cnt;
    int src, tag;
    int assigned_cnt = 0, terminated = 0;

    // Get filenames from directory
	vector<fs::path> filenames = get_names(argv[DIR_ARG], ".txt");
    file_cnt = filenames.size();
    
    // Buffer to hold profile vector from workers and 2d buffer to hold all vectors
    vector<char> buffer(300);
    vector<unordered_map<string, int>> profiles(file_cnt);
    
    // Keeps track of assigned files to each process
    vector<int> assigned(p);
    int profile_index = 0;

    do {
        // Get profile vector from worker
        MPI_Recv(&buffer[0], buffer.size(), MPI_BYTE, 
                MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        src = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        if(tag == VECTOR_MSG) {
            int size;
            MPI_Get_count(&status, MPI_BYTE, &size);    // Get buffer size 
            buffer.resize(size);    // Resize the buffer

            // The profile vector is serialized 
            // It is being deserialized here
            istringstream iss(string(&buffer[0], buffer.size()));
            archive::binary_iarchive ia(iss);
            unordered_map<string, int> umap;
            ia >> umap;
            profiles[profile_index++] = umap;
        }
        // Assign work or stop
        if(assigned_cnt < file_cnt) {
            const char * filename = filenames[assigned_cnt].c_str();
            MPI_Send(filename, strlen(filename) + 1, MPI_CHAR, src, FILE_NAME_MSG, MPI_COMM_WORLD);
            assigned[src] = assigned_cnt++;
        }
        // Send NULL to worker process, indicating termination
        else {
            MPI_Send(NULL, 0, MPI_CHAR, src, FILE_NAME_MSG, MPI_COMM_WORLD);
            terminated++;
        }
    } while(terminated < (p-1));

    // Printing the profile vectors
    int index = 0;
    for(auto inner_dict: profiles) {
        cout << "-----------------------------------" << endl;
        cout << filenames[index++] << endl;  
        for(auto j: inner_dict) {
                cout << j.first << ": " << j.second << " ";
            }
        cout << endl;
    }
}

void worker(int argc, char* argv[], MPI_Comm worker_comm) {
    MPI_Request pending;
    MPI_Status status;
    int worker_id;
    int dict_size, name_len;
    char* name;
    vector<char*> buffer;
    
    MPI_Comm_rank(worker_comm, &worker_id);
    // Send initial request for work
    MPI_Isend(NULL, 0, MPI_BYTE, 0, EMPTY_MSG, MPI_COMM_WORLD, &pending);

    // Read dictionary from file
    // TODO:: Optimize this by making process 0 read dictionary and broadcasting to other processes
    read_dictionary(argv[DICT_ARG], buffer, dict_size);

    while(1) {
        // Find length of file names
        MPI_Probe(0, FILE_NAME_MSG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &name_len);

        // Break out of loop if no work 
        if(!name_len) break;

		// Receive file name from manager
        name = (char*)malloc(name_len);
        MPI_Recv(name, name_len, MPI_CHAR, 0, FILE_NAME_MSG, MPI_COMM_WORLD, &status);

        // Creating and initializing hashtable dict, which would hold profile vector 
        unordered_map<string, int> dict;
        for(auto i: buffer) dict.insert({string(i), 0});

        // Create profile vector for filename
        make_profile(name, dict);
        free(name);

		// MPI does not support sending unordered maps hence it needs to be serialized
        // Serialize the dictionary
        ostringstream oss;
        archive::binary_oarchive oa(oss);
        oa << dict;

        // Send serialized profile vectors to manager
        MPI_Send(oss.str().c_str(), oss.str().size(), MPI_BYTE, 0, VECTOR_MSG, MPI_COMM_WORLD);
    }
}

void make_profile(char* filename, unordered_map<string, int>& dict)
{
    // Function to make profile vectors.
    fstream file;
    string word;

    file.open(filename);
    while(file >> word)
    {
        // if there is word in filename that is present in the dictionary then increment its count
        unordered_map<string, int>::iterator it = dict.find(word);
        if(it != dict.end())
            it->second++;
    }
    file.close();
}

void read_dictionary(char* filename, vector<char*>& buffer, int& file_len)
{
    // Function to read from dictionary file and write to buffer
    fstream file;
    string word;
    file_len = 0;

    file.open(filename);
    while(file >> word)
    {
        // convert string to char array before pushing into buffer
        char* cstr = new char[word.length()+1];
        strcpy(cstr, word.c_str());

        buffer.push_back(cstr);
        file_len++;
    }
}

vector<fs::path> get_names(fs::path const & root, string const & ext)
{
    // Function that returns path of all ext files in directory structure
    vector<fs::path> paths;

    if (fs::exists(root) && fs::is_directory(root))
    {
        for (auto const & entry : fs::recursive_directory_iterator(root))
        {
            if (fs::is_regular_file(entry) && entry.path().extension() == ext)
                paths.emplace_back(entry.path());
        }
    }
    return paths;
}        

#define BOOST_FILESYSTEM_VERSION 3
#define BOOST_FILESYSTEM_NO_DEPRECATED 
#include <boost/filesystem.hpp>
#include <bits/stdc++.h>
namespace fs = boost::filesystem;

/**
 *   Return the filenames of all files that have the specified extension
 *          in the specified directory and all subdirectories.

 *   https://stackoverflow.com/questions/11140483/how-to-get-list-of-files-with-a-specific-extension-in-a-given-folder
 */
std::vector<fs::path> get_all(fs::path const & root, std::string const & ext)
{
    std::vector<fs::path> paths;

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

int main() {
	std::vector<fs::path> paths = get_all("./", ".txt");
	for(auto i: paths) {
		std::cout << i.c_str() << std::endl;
	}
}
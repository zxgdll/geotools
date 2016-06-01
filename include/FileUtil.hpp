#ifndef __FILEUTIL_HPP__
#define __FILEUTIL_HPP__

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <vector>
#include <string>

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;

class FileUtil {
public:

	/**
	 * Populates the vector with the files contained in dir. If ext is specified, filters
	 * the files by that extension (case-insensitive). If dir is a file, it is added to the list.
	 * Returns the number of files found.
	 */
	static int dirlist(std::string &dir, std::vector<std::string> &files, std::string &ext = std::string()) {
		std::list<fs::path> files;
		if(fs::is_regular_file(dir)) {
			files.push_back(src);
		} else {
			fs::directory_iterator end;
			fs::directory_iterator di(src);
			for(; di != end; ++di) {
				if(!ext.empty()) {
					std::string p(di->path().string());
					alg::to_lower(p);
					if(alg::ends_with(p, ext))
						files.push_back(di->path());
				} else {
					files.push_back(di->path());
				}
			}
		}
		return files.size();
	}	

};

#endif
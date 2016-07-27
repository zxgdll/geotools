#ifndef __LASVALIDATE_HPP__
#define __LASVALIDATE_HPP__

#include <set>
#include <vector>
#include <string>

#include "geotools.h"

void validate(std::string &outfile, std::string &pointfile, std::string &datafile, std::vector<std::string> &lasfiles, 
	std::set<int> &classes, double distance);

#endif
#include "util.hpp"
#include <memory>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "liblas/liblas.hpp"

using namespace geotools::util;

int main(int argc, char **argv) {

	std::string lasfile = argv[1];

	/*
	std::list<double> lst;
	std::cerr << "size of list " << sizeof(lst) << std::endl;
	lst.push_back(131.0);
	std::cerr << "size of list with 1 " << sizeof(lst) << std::endl;
	lst.push_back(133.0);
	std::cerr << "size of list with 2 " << sizeof(lst) << std::endl;
	*/

	/*
	std::string file = Util::tmpFile();
	std::unique_ptr<MappedFile> mfile = Util::mapFile(file, sizeof(std::vector<double>) + sizeof(double) * 1000);
	void *ptr = mfile->data();
	std::vector<double> *lst = (std::vector<double>*) ptr;
	for(int i = 0; i < 1000; ++i)
		lst->push_back((double) (rand() % 1000));
	
	//std::sort(lst->begin(), lst->end());
	for(int i = 0; i < 1000; ++i)
		std::cerr << (*lst)[i] << std::endl;
	*/


	struct {
		bool operator()(const liblas::Point &p1, const liblas::Point &p2) {
			return p1.GetX() < p2.GetX();
		}
	} sortfn;

	liblas::ReaderFactory rf;
	std::ifstream instr(lasfile.c_str(), std::ios::in|std::ios::binary);
	liblas::Reader lasReader = rf.CreateWithStream(instr);
	liblas::Header lasHeader = lasReader.GetHeader();

	std::string file = Util::tmpFile();
	std::unique_ptr<MappedFile> mfile = Util::mapFile(file, 
		sizeof(std::vector<liblas::Point>) + sizeof(liblas::Point) * lasHeader.GetPointRecordsCount());
	void *ptr = mfile->data();
	std::vector<liblas::Point> *lst = (std::vector<liblas::Point>*) ptr;

	while(lasReader.ReadNextPoint()) {
		liblas::Point pt = lasReader.GetPoint();
		lst->push_back(std::move(pt));
	}

	std::sort(lst->begin(), lst->end(), sortfn);

	for(int i = 0; i < lst->size(); ++i)
		std::cerr << (*lst)[i].GetX() << ", " << (*lst)[i].GetY() << std::endl;
}
/**
 */

#include "geotools.h"
#include "Util.hpp"
#include "Raster.hpp"

#include "rasterfilter.hpp"

using namespace raster;

template <class T>
Link<T>::Link(T value) {
	this->value = value;
}

template <class T>
Chain<T>::Chain() :
	root(nullptr), current(nullptr),
	count(0) {
}

template <class T>
Link<T>* Chain<T>::get(int index) {
	if(index < 0 || index >= count)
		_argerr("Invalid index " << index);
	Link<T> *t = root;
	int i = 0;
	while(i++ < index && t->next)
		t = t->next;
	return t;
}

template <class T>
double Chain<T>::median() {
	if(count % 2 == 0) {
		Link<T> *a = get(count / 2 - 1);
		Link<T> *b = a->next;
		return (a->value + b->value) / 2.0;
	} else {
		return get(count / 2)->value;
	}
}

template <class T>
void Chain<T>::clear() {
	while(root) {
		current = root->next;
		delete root;
		root = current;
	}
	count = 0;
	current = root = nullptr;
}

template <class T>
Chain<T>::~Chain() {
	clear();
}

template <class T>
void Chain<T>::insert(T value) {
	if(!root) {
		current = root = new Link<T>(value);
		++count;
	} else {
		Link<T> *tmp = current;
		if(value < tmp->value) {
			while(true) {
				if(value < tmp->value) {
					if(tmp->prev) {
						tmp = tmp->prev;
					} else {
						root = tmp->prev = new Link<T>(value);
						tmp->prev->next = tmp;
						++count;
						break;
					}
				} else {
					Link<T> *n = new Link<T>(value);
					n->next = tmp;
					n->prev = tmp->prev;
					n->next->prev = n;
					n->prev->next = n;
					++count;
					break;
				}
			}
		} else {
			while(true) {
				if(value >= tmp->value) {
					if(tmp->next) {
						tmp = tmp->next;
					} else {
						tmp->next = new Link<T>(value);
						tmp->next->prev = tmp;
						++count;
						break;
					}
				} else {
					Link<T> *n = new Link<T>(value);
					n->prev = tmp;
					n->next = tmp->next;
					n->next->prev = n;
					n->prev->next = n;
					++count;
					break;
				}
			}
		}
	}
}

template <class T>
MedianFilter<T>::MedianFilter(int windowSize) :
	m_windowSize(windowSize) {
}

template <class T>
void MedianFilter<T>::process(Raster<T> &inrast, Raster<T> &outrast) {
	if(m_windowSize <= 0)
		_argerr("The window size must be larger than 1. " << m_windowSize << " given.");
	if(m_windowSize % 2 == 0) {
		_warn("Window size is even; adding 1.");
		++m_windowSize;
	}

	int count = m_windowSize * m_windowSize;
	int offset = m_windowSize / 2 + 1;
	Chain<T> chain;

	for(int r = offset; r < inrast.rows() - offset; ++r) {
		for(int c = offset; c < inrast.cols() - offset; ++c) {

			int i = 0;
			for(int rr = 0; rr < m_windowSize; ++rr) {
				for(int cc = 0; cc < m_windowSize; ++cc)
					chain.insert(inrast.get(cc, rr));
			}

			outrast.set(c, r, chain.median());
			chain.clear();

		}
	}
}

void usage() {
	std::cerr << "Usage: rasterfilter [options] <input file> <output file>\n"
		<< "    -f <filter>       -- The type of filter. Available types: median.\n"
		<< "    -w <window size>  -- For filters with a window size, gives the size in pixels.\n"
		<< "    -v               Verbose mode.\n";
}

int main(int argc, char **argv) {

 	try {

	 	int window = 0;
	 	std::vector<std::string> files;
	 	std::string filterName;
	 	bool verbose = false;

	 	for(int i = 1; i < argc; ++i) {
	 		std::string arg(argv[i]);
	 		if(arg == "-f") {
	 			filterName.assign(argv[++i]);
	 		} else if(arg == "-v") {
	 			verbose = true;
	 		} else if(arg == "-w") {
	 			window = atoi(argv[++i]);
	 		} else {
	 			files.push_back(argv[i]);
	 		}
	 	}

	 	if(files.size() < 2)
	 		_argerr("An input and output file are required.");

	 	std::string infile = files[0];
	 	std::string outfile = files[1];

	 	_loglevel(verbose ? LOG_TRACE : LOG_ERROR);

	 	Filter<float> *filter = nullptr;
	 	if(filterName == "median") {
	 		filter = new MedianFilter<float>(window);
	 	} else {
	 		_argerr("Unknown filter type: " << filterName);
	 	}

	 	Raster<float> inrast(infile);
	 	Raster<float> outrast(outfile, inrast);

	 	filter->process(inrast, outrast);

 	} catch(const std::exception &e) {
 		std::cerr << e.what() << std::endl;
 		usage();
 		return 1;
 	}

 	return 0;
 }

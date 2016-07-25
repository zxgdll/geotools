#ifndef __RASTERFILTER_HPP__
#define __RASTERFILTER_HPP__

#include <Raster.hpp>

using namespace raster;

template <class T>
class Link {
public:
	T value;
	Link<T> *next;
	Link<T> *prev;
	Link(T value);
};

template <class T>
class Chain {
public:
	Link<T> *root;
	Link<T> *current;
	int count;
	Chain();

	void insert(T value);

	Link<T>* get(int index);

	double median();

	void clear();

	~Chain();

};


template <class T>
class Filter {
public:
	virtual void process(Raster<T> &inrast, Raster<T> &outrast) =0;
};

template <class T>
class MedianFilter : public Filter<T> {
private:
	int m_windowSize;
public:
	MedianFilter(int windowSize);
	void process(Raster<T> &inrast, Raster<T> &outrast);

};

#endif
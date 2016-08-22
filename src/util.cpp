#include <set>
#include <list>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <map>

#include "csv.h"

#include "geotools.h"
#include "util.hpp"

using namespace geotools::util;

Point::Point(double x, double y, double z) :
	x(x), y(y), z(z) {
}

Point::Point(double x, double y, double z, const std::map<std::string, std::string> &fields) :
	Point(x, y, z) {
	for(auto it : fields)
		this->fields[it.first] = it.second;
}

Point::Point() : 
	Point(0, 0, 0) {
}

Bounds::Bounds() : Bounds(G_DBL_MAX_NEG, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_POS) {
}

Bounds::Bounds(double minx, double miny, double maxx, double maxy) :
	Bounds(minx, miny, maxx, maxy, G_DBL_MAX_NEG, G_DBL_MAX_POS) {
}

Bounds::Bounds(double minx, double miny, double maxx, double maxy, double minz, double maxz) :
	m_minx(minx), m_miny(miny), 
	m_maxx(maxx), m_maxy(maxy),
	m_minz(minz), m_maxz(maxz) {
}

bool Bounds::contains(double x, double y) const {
	return x >= m_minx && x <= m_maxx && y >= m_miny && y <= m_maxy;
}

bool Bounds::contains(double x, double y, double z) const {
	return contains(x, y) && z >= m_minz && z <= m_maxz;
}

bool Bounds::contains(const Bounds &b) const {
	return contains(b.minx(), b.miny(), b.minz()) && contains(b.maxx(), b.maxy(), b.maxz());
}

bool Bounds::intersects(const Bounds &b, int dims) const {
	if(dims == 3) {
		return b.contains(*this) || 
			contains(b.minx(), b.miny(), b.minz()) || contains(b.minx(), b.maxy(), b.minz()) || 
			contains(b.maxx(), b.miny(), b.minz()) || contains(b.maxx(), b.maxy(), b.minz()) ||
			contains(b.minx(), b.miny(), b.maxz()) || contains(b.minx(), b.maxy(), b.maxz()) || 
			contains(b.maxx(), b.miny(), b.maxz()) || contains(b.maxx(), b.maxy(), b.maxz());
	} else {
		return b.contains(*this) ||
			contains(b.minx(), b.miny()) || contains(b.minx(), b.maxy()) || 
			contains(b.maxx(), b.miny()) || contains(b.maxx(), b.maxy());
	}
}

double Bounds::minx() const {
	return m_minx;
}

void Bounds::minx(double minx) {
	m_minx = minx;
}

double Bounds::miny() const {
	return m_miny;
}

void Bounds::miny(double miny) {
	m_miny = miny;
}

double Bounds::minz() const {
	return m_minz;
}

void Bounds::minz(double minz) {
	m_minz = minz;
}

double Bounds::maxx() const {
	return m_maxx;
}

void Bounds::maxx(double maxx) {
	m_maxx = maxx;
}

double Bounds::maxy() const {
	return m_maxy;
}

void Bounds::maxy(double maxy) {
	m_maxy = maxy;
}

double Bounds::maxz() const {
	return m_maxz;
}

void Bounds::maxz(double maxz) {
	m_maxz = maxz;
}

double Bounds::width() const {
	return maxx() - minx();
}

double Bounds::height() const {
	return maxy() - miny();
}

double Bounds::depth() const {
	return maxz() - minz();
}

int Bounds::cols(double resolution) const {
	return g_max(1, (int) std::ceil(width() / resolution));
}

int Bounds::rows(double resolution) const {
	return g_max(1, (int) std::ceil(height() / resolution));
}

void Bounds::extend(const Bounds &b) {
	m_minx = g_min(b.minx(), m_minx);
	m_maxx = g_max(b.maxx(), m_maxx);
	m_miny = g_min(b.miny(), m_miny);
	m_maxy = g_max(b.maxy(), m_maxy);
	m_minz = g_min(b.minz(), m_minz);
	m_maxz = g_max(b.maxz(), m_maxz);
}

void Bounds::extendX(double x) {
	m_minx = g_min(x, m_minx);
	m_maxx = g_max(x, m_maxx);
}

void Bounds::extendY(double y) {
	m_miny = g_min(y, m_miny);
	m_maxy = g_max(y, m_maxy);
}

void Bounds::extendZ(double z) {
	m_minz = g_min(z, m_minz);
	m_maxz = g_max(z, m_maxz);
}

void Bounds::extend(double x, double y) {
	extendX(x);
	extendY(y);
}

void Bounds::extend(double x, double y, double z) {
	extend(x, y);
	extendZ(z);
}

double Bounds::operator[](size_t pos) const {
	switch(pos) {
	case 0: return m_minx;
	case 1: return m_miny;
	case 2: return m_maxx;
	case 3: return m_maxy;
	case 4: return m_minz;
	case 5: return m_maxz;
	default:
		g_argerr("Illegal position: " << pos);
	}
}

void Bounds::snap(double resolution) {
	minx(std::floor(minx() / resolution) * resolution);
	miny(std::floor(miny() / resolution) * resolution);
	maxx(std::floor(maxx() / resolution) * resolution + resolution);
	maxy(std::floor(maxy() / resolution) * resolution + resolution);
}

void Bounds::collapse(int dims) {
	minx(G_DBL_MAX_POS);
	miny(G_DBL_MAX_POS);
	maxx(G_DBL_MAX_NEG);
	maxy(G_DBL_MAX_NEG);
	if(dims == 3) {
		minz(G_DBL_MAX_POS);
		maxz(G_DBL_MAX_NEG);
	}
}
	
std::string Bounds::print() const {
	std::stringstream s;
	print(s);
	return s.str();
}

void Bounds::print(std::ostream &str) const {
	str << "[Bounds: " << minx() << ", " << miny() << ", " << minz() << "; " << maxx() << ", " << maxy() << ", " << maxz() << "]";
}

/**
 * Split a comma-delimited string into a set of unique integers.
 */
void Util::intSplit(std::set<int> &values, const char *str) {
	std::stringstream ss(str);
	std::string item;
	while(std::getline(ss, item, ','))
		values.insert(atoi(item.c_str()));
}

/**
 * Split a comma-delimited string into a set of unique integers.
 */
void Util::intSplit(std::list<int> &values, const char *val) {
	std::stringstream ss(val);
	std::string item;
	while(std::getline(ss, item, ','))
		values.push_back(atoi(item.c_str()));
}

/**
 * Split a comma-delimited string into a set of unique integers.
 */
void Util::intSplit(std::vector<int> &values, const char *str) {
	std::stringstream ss(str);
	std::string item;
	while(std::getline(ss, item, ','))
		values.push_back(atoi(item.c_str()));
}

/**
 * Return true if the integer is in the set, or the set is empty.
 */
bool Util::inList(std::set<int> &values, int value) {
	return values.size() == 0 || values.find(value) != values.end();
}

bool Util::inList(std::vector<int> &values, int value) {
	return std::find(values.begin(), values.end(), value) != values.end();
}

void Util::copyfile(std::string &srcfile, std::string &dstfile) {
	std::ifstream src(srcfile.c_str(), std::ios::binary);
	std::ofstream dst(dstfile.c_str(), std::ios::binary);
	dst << src.rdbuf();
}

/**
 * Load the samples from a csv file. The file must have x, y and z headers.
 */
void Util::loadXYZSamples(std::string &datafile, std::vector<std::tuple<double, double, double> > &samples) {
	io::CSVReader<3> in(datafile.c_str());
	in.read_header(io::ignore_extra_column, "x", "y", "z");
	double x, y, z;
	while(in.read_row(x, y, z))
		samples.push_back(std::make_tuple(x, y, z));
}

void Util::loadIDXYZSamples(std::string &datafile, std::vector<std::tuple<std::string, double, double, double> > &samples) {
	io::CSVReader<4> in(datafile.c_str());
	in.read_header(io::ignore_extra_column, "id", "x", "y", "z");
	std::string id;
	double x, y, z;
	while(in.read_row(id, x, y, z))
		samples.push_back(std::make_tuple(id, x, y, z));
}

void Util::status(int step, int of, const std::string &message, bool end) {
        #pragma omp critical(__status)
        {
                if(step < 0)  step = 0;
                if(of <= 0)   of = 1;
                if(step > of) of = step;
                float status = (float) (step * 100) / of;
                std::stringstream out;
                out << "Status: " << std::fixed << std::setprecision(2) << status << "% " << message << std::right << std::setw(100) << std::setfill(' ');
                if(end) 
			out << std::endl;
		else
			out << '\r';
		std::cerr << out.str();
                std::cerr.flush();
        }
}


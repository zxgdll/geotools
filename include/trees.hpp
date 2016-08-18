#ifndef __TREETOPS_HPP__
#define __TREETOPS_HPP__

#include <string>
#include <map>

namespace trees {

	namespace util {

	/**
                 * A simple class for maintaining information about a tree top.
                 */
                class Top {
                public:
                        size_t m_id;
                        double m_x;
                        double m_y;
                        double m_z;
                        int m_col;
                        int m_row;
                        Top(size_t id, double x, double y, double z, int col, int row) :
                                m_id(id),
                                m_x(x), m_y(y), m_z(z),
                                m_col(col), m_row(row) {
                        }
                        Top(const Top &t) :
                                Top(t.m_id, t.m_x, t.m_y, t.m_z, t.m_col, t.m_row) {
                        }
                        Top() :
                                Top(0, 0, 0, 0, 0, 0) {
                        }
                };

	}
	
	void treetops(const std::string &inraster, const std::string &topsvect, std::map<size_t, std::unique_ptr<trees::util::Top> > &tops,  
		int window, double minHeight, const std::string &smoothed, double sigma = 0.8408964, int size = 3);

 	void treecrowns(const std::string &infile, const std::string &outrfile, const std::string &outvfile,
                std::map<size_t, std::unique_ptr<trees::util::Top> > &tops, double threshold, double radius);
}

#endif

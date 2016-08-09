#include <set>
#include <list>
#include <sstream>
#include <algorithm>
#include <fstream>

#include <liblas/liblas.hpp>

#include "csv.h"

#include "geotools.h"

namespace geotools {

	namespace util {

		class Bounds {
		private:
			double m_minx, m_miny, m_minz;
			double m_maxx, m_maxy, m_maxz;
		public:
			Bounds();

			Bounds(double minx, double miny, double maxx, double maxy);

			Bounds(double minx, double miny, double maxx, double maxy, double minz, double maxz);

			bool contains(double x, double y);

			bool contains(double x, double y, double z);

			bool contains(geotools::util::Bounds &b);

			bool intersects(Bounds &b, int dims = 2);

			double minx();

			void minx(double minx);

			double miny();

			void miny(double miny);

			double minz();

			void minz(double minz);

			double maxx();

			void maxx(double maxx);

			double maxy();

			void maxy(double maxy);

			double maxz();

			void maxz(double maxz);

			double width();

			double height();

			double depth();

			int cols(double resolution);

			int rows(double resolution);

			void extend(Bounds &b);

			void extendX(double x);
			
			void extendY(double y);
			
			void extendZ(double z);
			
			void extend(double x, double y);
			
			void extend(double x, double y, double z);

			void collapse(int dims = 2);

			double operator[](size_t pos);

			void snap(double resolution);
	
			std::string print();

			void print(std::ostream &str);
		
		};

		/**
		 * Provides utility methods for working with LiDAR data.
		 */
		class Util {
		public:

		#ifdef LIBLAS_HPP_INCLUDED

			/**
			 * Compute the bounds of a LAS file using the Header and a double
			 * array to contain the result. If the dims argument is 2, the
			 * horizontal coords are considered and if it's 3, all coords are.
			 */
			static bool computeLasBounds(liblas::Header &hdr, geotools::util::Bounds &bounds, int dims = 2);

			/**
			 * Compute the bounds of a LAS file using the Reader and a double
			 * array to contain the result. If the dims argument is 2, the
			 * horizontal coords are considered and if it's 3, all coords are.
			 * Use this method when the header bounds are bogus (it iterates over all the points.
			 */
			static bool computeLasBounds(liblas::Reader &rdr, geotools::util::Bounds &bounds, int dims = 2);

		#endif

			
			/**
			 * Split a comma-delimited string into a set of unique integers.
			 */
			static void intSplit(std::set<int> &values, const char *str);
		

			/**
			 * Split a comma-delimited string into a set of unique integers.
			 */
			static void intSplit(std::list<int> &values, const char *val);

			/**
			 * Split a comma-delimited string into a set of unique integers.
			 */
			static void intSplit(std::vector<int> &values, const char *str);

			/**
			 * Return true if the integer is in the set, or the set is empty.
			 */
			static bool inList(std::set<int> &values, int value);

			static bool inList(std::vector<int> &values, int value);

			/**
			 * Prints out a status message; a percentage representing current
			 * of total steps.
			 */
			static void status(int current, int total);

			static void copyfile(std::string &srcfile, std::string &dstfile);

			/**
			 * Load the samples from a csv file. The file must have x, y and z headers.
			 */
			static void loadXYZSamples(std::string &datafile, std::vector<std::tuple<double, double, double> > &samples);

			static void loadIDXYZSamples(std::string &datafile, std::vector<std::tuple<std::string, double, double, double> > &samples);

			static void status(int step, int of, bool end = false);

		};

	} // util

} // geotools

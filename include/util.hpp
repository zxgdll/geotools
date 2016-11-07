#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <set>
#include <list>
#include <ostream>
#include <vector>
#include <map>
#include <memory>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>

namespace geotools {

	namespace util {


		// Provides access to an allocated buffer
		// which will be safely disposed of.
		class Buffer {
		public:
			void *buf;
			Buffer(uint64_t size) {
				buf = std::calloc(size, 1);
			}
			~Buffer() {
				std::free(buf);
			}
		};

		class Callbacks {
		public:
			virtual ~Callbacks()=0;
			virtual void fileCallback(float status) const =0;
			virtual void overallCallback(float status) const =0;
		};

		class Point {
                public:
			double x, y, z;
                        std::map<std::string, std::string> fields;
			Point();
			Point(double x, double y, double z = 0);
			Point(double x, double y, double z, const std::map<std::string, std::string> &fields);
                };

		class Bounds {
		private:
			double m_minx, m_miny, m_minz;
			double m_maxx, m_maxy, m_maxz;
		public:
			Bounds();

			Bounds(double minx, double miny, double maxx, double maxy);

			Bounds(double minx, double miny, double maxx, double maxy, double minz, double maxz);

			bool contains(double x, double y) const;

			bool contains(double x, double y, double z) const;

			bool contains(const geotools::util::Bounds &b) const;

			bool intersects(const geotools::util::Bounds &b, int dims = 2) const;

			Bounds intersection(const Bounds &other) const;

			double minx() const;

			void minx(double minx);

			double miny() const;

			void miny(double miny);

			double minz() const;

			void minz(double minz);

			double maxx() const;

			void maxx(double maxx);

			double maxy() const;

			void maxy(double maxy);

			double maxz() const;

			void maxz(double maxz);

			double width() const;

			double height() const;

			double depth() const;

			int cols(double resolution) const;

			int rows(double resolution) const;

			void extend(const geotools::util::Bounds &b);

			void extendX(double x);
			
			void extendY(double y);
			
			void extendZ(double z);
			
			void extend(double x, double y);
			
			void extend(double x, double y, double z);

			void collapse(int dims = 2);

			double operator[](size_t pos) const;

			void snap(double resolution);
	
			std::string print() const;

			void print(std::ostream &str) const;
		
		};

		class Util;

		/**
		 * Maintains a memory-mapped file, and gives access to the mapped data.
		 */
		class MappedFile {
			friend class Util;
		private:
			std::string m_filename;
			uint64_t m_size;
			boost::interprocess::file_mapping *m_mapping;
			boost::interprocess::mapped_region *m_region;
			bool m_remove;
		protected:
			MappedFile(const std::string &filename, uint64_t size, bool remove);
		public:
			void* data();
			uint64_t size();
			~MappedFile();
		};


		/**
		 * Provides utility methods for working with LiDAR data.
		 */
		class Util {
		public:
			
			static void parseRanges(std::set<double> &values, const char *str, double step = 1.0);

			static void parseRanges(std::set<int> &values, const char *str);

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

			static void intSplit(std::set<uint8_t> &values, const char *str);
			/**
			 * Return true if the integer is in the set, or the set is empty.
			 */
			static bool inList(std::set<int> &values, int value);

			static bool inList(std::vector<int> &values, int value);

			static void splitString(const std::string &str, std::list<std::string> &lst);
			// TODO: Use back inserter.
			static void splitString(const std::string &str, std::vector<std::string> &lst);
			
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

			static void status(int step, int of, const std::string &message = "", bool end = false);

			static const std::string tmpFile(const std::string &root);
			static const std::string tmpFile();
			static bool rm(const std::string &name);
			static bool mkdir(const std::string &dir);

			static std::unique_ptr<MappedFile> mapFile(const std::string &filename, uint64_t size, bool remove = true);

		};

	} // util

} // geotools

#endif

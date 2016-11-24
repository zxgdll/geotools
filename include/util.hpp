#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <set>
#include <list>
#include <fstream>
#include <vector>
#include <map>
#include <memory>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <cmath>
#include <unordered_map>
#include <sstream>
#include <string>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>

#include "geotools.h"

#ifdef _MSC_VER
#include <float.h>
namespace std {

    inline bool isnan(double value) {
        return _isnan(value);
    }
    
    inline double round(double value) {
        return value < 0 ? -std::floor(0.5 - value) : std::floor(0.5 + value);
    }
}
#endif

namespace geotools {

    namespace util {

        template <class T>
        class blocking_queue {
        private:
            std::queue<T> m_q;
            std::condition_variable m_c;
            std::mutex m_m;
        public:

            bool pop(T *idx) {
                bool ret = false;
                {
                    std::unique_lock<std::mutex> lk(m_m);
                    m_c.wait(lk);
                    if (!m_q.empty()) {
                        *idx = m_q.front();
                        m_q.pop();
                        ret = true;
                    }
                }
                return ret;
            }

            void push(T idx) {
                {
                    std::lock_guard<std::mutex> lk(m_m);
                    m_q.push(idx);
                }
                m_c.notify_one();
            }

            size_t size() {
                std::lock_guard<std::mutex> lk(m_m);
                return m_q.size();
            }

            bool empty() {
                std::lock_guard<std::mutex> lk(m_m);
                return m_q.empty();
            }

            void flush() {
                while (!empty())
                    m_c.notify_one();
            }

            void finish() {
                m_c.notify_all();
            }
        };

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

        // Provides methods for handling status callbacks.
        class Callbacks {
        public:
            virtual ~Callbacks() = 0;
            virtual void stepCallback(float status) const = 0;
            virtual void overallCallback(float status) const = 0;
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

            // Populates the vector with the files contained in dir. If ext is specified, filters
            // the files by that extension (case-insensitive). If dir is a file, it is added to the list.
            // Returns the number of files found.
            static int dirlist(const std::string &dir, std::vector<std::string> &files, 
                const std::string &ext = std::string());

            static std::unique_ptr<MappedFile> mapFile(const std::string &filename, 
                uint64_t size, bool remove = true);

        };

    } // util

} // geotools

#endif

#ifndef __SORTEDPOINTSTREAM_HPP__
#define __SORTEDPOINTSTREAM_HPP__

#include <map>
#include <string>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <atomic>
#include <queue>

#include "liblas/liblas.hpp"

#include "geos/geom/Geometry.h"

#include "util.hpp"

using namespace geotools::util;

namespace geotools {

    namespace las {

        class LASPoint {
        private:

            // Determine the scale factors used by all
            // points (coordinates are stored as ints in 
            // the LAS format.)
            //static double m_scaleX, m_scaleY, m_scaleZ;

            // Indicates the size taken up by a single
            // LASPoint in memory or on disc.
            const static uint64_t m_dataSize = 20;

        public:

            double x, y, z;
            uint16_t intensity;
            uint16_t returnNum, numReturns;
            uint16_t scanDirection;
            uint8_t cls;
            int8_t scanAngle;

            LASPoint();

            LASPoint(const liblas::Point &pt);

            ~LASPoint();

            // Sets the scale values used by all points.
            static void setScale(double x, double y, double z);

            // Set fields from a liblas Point object.
            void update(const liblas::Point &pt);

            bool operator<(const LASPoint&) const;
            bool operator>(const LASPoint&) const;
            bool operator==(const LASPoint&) const;
            bool operator!=(const LASPoint&) const;
            bool operator<=(const LASPoint&) const;
            bool operator>=(const LASPoint&) const;

            void write(std::ostream &str) const;

            void read(std::istream &str);

            void write(std::FILE *str);

            void read(std::FILE *str);

            void write(void *str) const;

            void read(void *str);

            // Return true if this is a last return.
            bool last() const;

            // Return true if this is a first return.
            bool first() const;

            // Return true if this is an intermediate return.
            bool intermediate() const;

            // Return true if this is a ground point (class 2).
            bool ground() const;

            // Return true if this is a single return.
            bool single() const;

            static double scaleX();
            static double scaleY();
            static double scaleZ();

            static uint64_t dataSize();

        };

        // Provides access to points from a set of
        // LAS (or other) files which have been ordered
        // into rows and cached as a file. This is not a
        // complete ordering, just one that makes further
        // processing easier.

        class SortedPointStream {
        private:
            // The list of source files.
            std::list<std::string> m_files;

            // The name of the file containing sorted points.
            std::string m_cacheFile;

            // The height of a row in map unites. May be negative
            // if the projection required it (see GDAL's affine transform.)
            double m_blockSize;

            // If true, rebuilds the cache file; otherwise reuses it
            // if it exists.
            bool m_rebuild;

            // Set to true after initialization.
            bool m_inited;

            // The file handle used for reading and writing.
            // Stays open between calls to next, but closes
            // on destruction.
            std::FILE *m_file;

            // If true, the working boundary is snapped to
            // the resolution.
            bool m_snap;

            // The number of threads to use for sorting
            // points.
            uint32_t m_threads;

            // The working bounds of the fileset.
            Bounds m_bounds;

            // The number of points in the fileset.
            uint64_t m_pointCount;

            // The number of rows in the cache.
            uint64_t m_rowCount;

            // The current row; used by next().
            uint64_t m_row;

            // Cache -- stores rows of points for writing; indexed by row.
            std::unordered_map<uint64_t, std::list<std::shared_ptr<LASPoint> > > m_cache;

            // Jump table -- stores the next offset to use for points 
            // associated with a given row.
            std::unordered_map<uint64_t, uint32_t> m_jump;

            // True if the consumer should be running.
            bool m_running;

            // The next available offset for writing rows.
            uint64_t m_nextJump;

            std::mutex m_fmtx;
            std::mutex m_cmtx;
            std::mutex m_wmtx;

            // The queue of files to consume.
            std::queue<std::string> m_fileq;

            // The queue of rows that need to be written
            // to disc.
            std::queue<uint64_t> m_flush;

            std::unique_ptr<geotools::util::MappedFile> m_mmap;

            LASPoint *m_lst;
            std::unique_ptr<MappedFile> m_mfile;
            uint64_t m_idx;

            // Read source files and populate the point cache.
            // To be run using (a) thread(s).
            void produce();

            // Write point cache to disc.
            // To be run using (a) thread(s).
            void consume();

        public:

            /**
            // Constructs a new point stream on the given files.
             *
            // The block size gives the size, in map units, of a block,
            // blocks are square and stored as individual files.
             */
            SortedPointStream(const std::list<std::string> &files, const std::string &cacheFile,
                    double blockSize, bool rebuild = true, bool snap = false, uint32_t threads = 1);

            // Initialize the point stream. This is a time-consuming 
            // operation.
            void init();

            // Set a callbacks object to handle progress.
            void setCallbacks(Callbacks *callbacks);

            ~SortedPointStream();

            uint64_t pointCount() const;

            uint64_t rowCount() const;

            uint64_t bufferSize() const;

            const Bounds& bounds() const;

            // Returns the next available row of points.
            bool next(std::list<std::shared_ptr<LASPoint> > &pts);
        };

    }
}


#endif

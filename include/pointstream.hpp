#ifndef __POINTSTREAM_HPP__
#define __POINTSTREAM_HPP__

#include <string>
#include <fstream>

#include <liblas/liblas.hpp>

#include "util.hpp"

namespace geotools {

    namespace las {

        class LASPoint {
        public:
            static double scaleX, scaleY, scaleZ;
            double x, y, z;
            uint16_t intensity;
            uint16_t ret, numRets;
            uint16_t scanDir;
            uint8_t cls;
            int8_t angle;

            LASPoint() :
            ret(0), numRets(0) {
            }

            LASPoint(const liblas::Point &pt) : LASPoint() {
                update(pt);
            }

            static void setScale(double x, double y, double z) {
                scaleX = x;
                scaleY = y;
                scaleZ = z;
            }

            void update(const liblas::Point &pt) {
                x = pt.GetX();
                y = pt.GetY();
                z = pt.GetZ();
                intensity = pt.GetIntensity();
                ret = pt.GetReturnNumber();
                numRets = pt.GetNumberOfReturns();
                cls = pt.GetClassification().GetClass();
                angle = pt.GetScanAngleRank();
            }

            bool last() const {
                return numRets > 0 && ret == numRets;
            }

            bool first() const {
                return numRets > 0 && ret == 1;
            }

            bool intermediate() const {
                return numRets > 2 && ret > 1 && ret < numRets;
            }

            bool ground() const {
                return cls == 2;
            }

            bool single() const {
                return numRets == 1;
            }
        };

        class PointStream {
        private:
            unsigned int m_fileIdx;
            unsigned long m_pointCount;
            unsigned long m_curPt;
            unsigned long m_curPtCount;
            geotools::util::Bounds m_bounds;
            std::unordered_map<std::string, geotools::util::Bounds> m_fileBounds;
            std::list<geotools::util::Bounds> m_completedBounds;
            std::vector<std::string> m_files;
            std::string m_currentFile;
            liblas::Reader *m_lasReader;
            std::ifstream *m_instr;

            /**
             * Initializes the point stream.
             */
            void init(const std::list<std::string> &files, bool deepBounds);

            bool loadNextReader();

        public:

            // TODO: Learn to use back_inserter
            static geotools::util::Bounds filterFiles(const std::list<std::string> &infiles, std::vector<std::string> &outfiles, const geotools::util::Bounds &bounds);

            /**
             * Constructs a new point stream on the given files.
             * Deep bounds determines whether the bounds of a file are computed from
             * the header alone, or from the points, should the header be incorrect.
             */
            PointStream(const std::list<std::string> &files, bool deepBounds = false);

            ~PointStream();

            geotools::util::Bounds bounds();

            unsigned int pointCount();

            bool contains(double x, double y);

            bool contains(double x1, double y1, double x2, double y2);

            bool containsCompleted(double x, double y);

            bool containsCompleted(double x1, double y1, double x2, double y2);
            /**
             * Reads the next available point into the LASPoint object.
             */
            bool next(LASPoint &pt, bool *lastPoint);

        };

    } // las

} // geotools


#endif
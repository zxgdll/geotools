#ifndef __LASPOINT_HPP__
#define __LASPOINT_HPP__

#include <fstream>
#include <cstdio>

#include "liblas/liblas.hpp"

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

            void readLAS0(char *buf);
            void readLAS1(char *buf);
            void readLAS2(char *buf);
            void readLAS6(char *buf);
            void readLAS8(char *buf);

        public:

            double x, y, z;
            uint16_t intensity;
            uint16_t returnNum, numReturns;
            uint16_t scanDirection;
            uint8_t cls;
            uint8_t clsFlags;
            int8_t scanAngle;
            uint8_t scanDir;
            bool isEdge;
            uint8_t red, green, blue, nir;
            uint8_t channel;
            double gpsTime;
            uint16_t sourceId;
            uint16_t userData;

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

            void readLAS(char *buf, uint8_t format);

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
    }
}


#endif

#include "lasutil.hpp"

using namespace geotools::util;

/**
 * Compute the bounds of a LAS file using the Header and a double
 * array to contain the result. If the dims argument is 2, the
 * horizontal coords are considered and if it's 3, all coords are.
 */
bool LasUtil::computeLasBounds(liblas::Header &hdr, geotools::util::Bounds &bounds, int dims) {
        if((dims != 2 && dims != 3) || hdr.GetMaxX() == hdr.GetMinX() || hdr.GetMaxY() == hdr.GetMinY())
                return false;
        bounds.collapse(dims == 3);
        if(dims == 2) {
                bounds.extend(hdr.GetMinX(), hdr.GetMinY());
                bounds.extend(hdr.GetMaxX(), hdr.GetMaxY());
        } else {
                bounds.extend(hdr.GetMinX(), hdr.GetMinY(), hdr.GetMinZ());
                bounds.extend(hdr.GetMaxX(), hdr.GetMaxY(), hdr.GetMaxZ());
        }
        return true;
}

/**
 * Compute the bounds of a LAS file using the Reader and a double
 * array to contain the result. If the dims argument is 2, the
 * horizontal coords are considered and if it's 3, all coords are.
 * Use this method when the header bounds are bogus (it iterates over all the points.
 */
bool LasUtil::computeLasBounds(liblas::Reader &rdr, geotools::util::Bounds &bounds, int dims) {
        if(dims != 2 && dims != 3)
                return false;
        bounds.collapse(dims == 3);
        while(rdr.ReadNextPoint()) {
                liblas::Point pt = rdr.GetPoint();
                if(dims == 2) {
                        bounds.extend(pt.GetX(), pt.GetY());
                } else {
                        bounds.extend(pt.GetX(), pt.GetY(), pt.GetZ());
                }
        }
        return true;
}



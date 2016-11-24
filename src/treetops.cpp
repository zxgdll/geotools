/*
 * treetops
 * 
 * This library provides methods for isolating tree tops and crowns from a
 * LiDAR (or other) canopy height model (CHM.) The output from this program is 
 * ideal for use with the spectral extraction module (spectral) or other 
 * analysis.
 * 
 * The usual sequence for producing useful output is:
 * 1) Smooth the original CHM. The CHM should possibly have pit- or 
 *    spike-removal applied to it first. Smoothing is a Gaussian kernel with
 *    configurable sigma and window size.
 * 2) Locate treetops. This is performed on the smoothed CHM and uses a 
 *    maximum-value kernel to locate local maxima. Locates the tree top height
 *    from the original CHM.
 * 3) Delineate crowns. Uses the tree tops as seeds in a region-growing 
 *    algorithm that creates tree crown boundaries with configurable limits.
 *
 *  Created on: May 3, 2016
 *      Author: Rob Skelly 
 *       Email: rob@dijital.ca
 */

#include <queue>
#include <iostream>
#include <atomic>

#include <omp.h>

#include "sqlite.hpp"
#include "util.hpp"
#include "raster.hpp"
#include "treetops.hpp"

using namespace geotools::raster;
using namespace geotools::util;
using namespace geotools::db;

using namespace geotools::treetops::config;
using namespace geotools::treetops::util;
using namespace geotools::treetops;

namespace geotools {

    namespace treetops {

        namespace util {

            // Represents a grid cell, and maintains some properties of the 
            // seed that originated it.
            class Node {
            public:
                uint64_t id;
                int c, r, tc, tr;
                double z, tz;

                Node(uint64_t id, int c, int r, double z, int tc, int tr, double tz) :
                    id(id),
                    c(c), r(r), tc(tc), tr(tr),
                    z(z), tz(tz) {
                }

                Node(const Top &top) :
                    id(top.id),
                    c(top.col), r(top.row), tc(top.col), tr(top.row),
                    z(top.z), tz(top.z) {
                }
            };

            // Returns true if the pixel at the center of the given raster is 
            // the maximum value in the raster.
            bool isMaxCenter(MemRaster<float> &raster, int col, int row, int window, double *max) {
                int cc = col + window / 2;
                int cr = row + window / 2;
                float nd = raster.nodata();
                if (raster.get(cc, cr) == nd)
                    return false;
                *max = 0;
                int mc = 0, mr = 0;
                for (int r = row; r < row + window; ++r) {
                    for (int c = col; c < col + window; ++c) {
                        float v = raster.get(c, r);
                        if (v != nd && v > *max) {
                            *max = v;
                            mc = c;
                            mr = r;
                        }
                    }
                }
                return mc == cc && mr == cr;
            }

            // Guess a value for a cell, based on its neighbours.
            double interpNodata(Grid<float> &rast, int col, int row) {
                int size = 1;
                double nodata = rast.nodata();
                double v, t;
                int n;
                while (size < 1000) {
                    n = 0;
                    t = 0;
                    for (int c = g_max(0, col - size); c < g_min(rast.cols(), col + size + 1); ++c) {
                        v = rast.get(c, g_max(0, row - size));
                        if (v != nodata) t += v, ++n;
                        v = rast.get(c, g_min(rast.rows() - 1, row + size));
                        if (v != nodata) t += v, ++n;
                    }
                    for (int r = g_max(1, row - size); r < g_min(rast.rows(), row + size + 1); ++r) {
                        v = rast.get(g_max(0, col - size), r);
                        if (v != nodata) t += v, ++n;
                        v = rast.get(g_min(rast.cols() - 1, col + size), r);
                        if (v != nodata) t += v, ++n;
                    }
                    if (n > 0)
                        return t / n;
                    ++size;
                }
                g_runerr("Couldn't find a pixel to use as fill.");
            }

        } // util

    } // trees

} // geotools

Top::Top(uint64_t id, double x, double y, double z, double uz, int col, int row) :
    id(id),
    x(x), y(y), z(z), uz(uz),
    col(col), row(row) {
}

Top::Top() :
    id(0),
    x(0), y(0), z(0), uz(0),
    col(0), row(0) {
}

void Treetops::setCallbacks(Callbacks *callbacks) {
    m_callbacks = callbacks;
}

void Treetops::smooth(const TreetopsConfig &config) {
    config.checkSmoothing();
    Raster<float> in(config.smoothOriginalCHM);
    Raster<float> out(config.smoothSmoothedCHM, 1, in);
    in.smooth(out, config.smoothStdDev, config.smoothWindowSize, m_callbacks);
}

void Treetops::treetops(const TreetopsConfig &config) {
    config.checkTops();
    if(m_callbacks)
        m_callbacks->stepCallback(0.01);

    // Initialize input rasters.
    g_debug(" -- tops: opening rasters");
    Raster<float> original(config.topsOriginalCHM);
    Raster<float> smoothed(config.topsSmoothedCHM);

    // Prepare database.
    g_debug(" -- tops: preparing db");
    std::map<std::string, int> fields;
    fields["id"] = 1;
    SQLite db(config.topsTreetopsDatabase, SQLite::POINT, config.srid, fields, true);
    db.makeFast();
    db.dropGeomIndex();
    db.setCacheSize(config.tableCacheSize);

    if(m_callbacks)
        m_callbacks->stepCallback(0.02);

    int32_t bufSize = 256;              // The number of rows each thread works on at a time.
    std::atomic<uint64_t> topCount(0);
    std::atomic<int32_t> curRow(0);     // Used for status indicator.
    std::atomic<uint64_t> topId(0);

    #pragma omp parallel
    {
        g_debug(" -- tops: processing " << omp_get_thread_num());
        MemRaster<float> blk(original.cols(), bufSize + config.topsWindowSize);
        blk.nodata(original.nodata());
        blk.fill(blk.nodata());

        std::map<uint64_t, std::unique_ptr<Top> > tops0;
        int32_t bufSize0;

        #pragma omp for
        for (int32_t brow = 0; brow < original.rows(); brow += bufSize) {

            bufSize0 = g_min(bufSize, original.rows() - brow - config.topsWindowSize);
            if(bufSize0 < bufSize)
                blk.init(original.cols(), bufSize0 + config.topsWindowSize);

            g_debug(" -- bufsize " << bufSize0);

            #pragma omp critical(__c)
            smoothed.readBlock(0, brow, blk);

            for (int32_t row = 0; row < bufSize0; ++row) {
                for (int32_t col = 0; col < original.cols() - config.topsWindowSize; ++col) {
        
                    int32_t r = row + config.topsWindowSize / 2;
                    int32_t c = col + config.topsWindowSize / 2;            
                    double max;

                    if (blk.get(c, r) >= config.topsMinHeight &&
                            isMaxCenter(blk, col, row, config.topsWindowSize, &max)) {

                        // Compute the id based on the cell.
                        uint64_t id = ((uint64_t) c << 32) | r;
                        // Get the original height from the unsmoothed raster.
                        double umax = original.get(c, r + brow);

                        std::unique_ptr<Top> pt(new Top(
                                ++topId,
                                original.toCentroidX(c), // center of pixel
                                original.toCentroidY(r + brow),
                                max,
                                umax,
                                c,
                                r + brow
                        ));
                        tops0[id] = std::move(pt);
                    }
                }
            }

            if(m_callbacks) {
                curRow += bufSize0;
                m_callbacks->stepCallback(0.02 + (float) curRow / original.rows() * 0.48);
            }   
        }

        uint64_t topCount0 = tops0.size();
        uint64_t b = 0;
        uint64_t batch = db.maxAddPointCount();

        topCount += topCount0;

        std::vector<std::unique_ptr<Point> > points;
        for (const auto &it : tops0) {
            const uint64_t &id = it.first;
            const std::unique_ptr<Top> &t = it.second;
            std::map<std::string, std::string> fields;
            fields["id"] = std::to_string(id);
            std::unique_ptr<Point> pt(new Point(t->x, t->y, t->uz, fields));
            points.push_back(std::move(pt));

            if (++b % batch == 0 || b >= topCount0) {
                g_debug("inserting " << points.size() << " points");
                #pragma omp critical(__b)
                {
                    db.begin();
                    db.addPoints(points);
                    db.commit();
                }
                points.clear();
                if(m_callbacks)
                    m_callbacks->stepCallback((float) b / topCount0 * 0.3 + 0.48);
            }
        }
    }

    if(m_callbacks)
        m_callbacks->stepCallback(0.99);

    if (config.buildIndex) {
        g_debug(" -- tops: building index");
        db.createGeomIndex();
    }
    db.makeSlow();

    g_debug(" -- tops: done.");
    if(m_callbacks)
        m_callbacks->stepCallback(1.0);

}

void Treetops::treecrowns(const TreetopsConfig &config) {
    config.checkCrowns();
    if(m_callbacks)
        m_callbacks->stepCallback(0.01);

    g_debug(" -- crowns: delineating crowns...");
    uint32_t batchSize = 10000;

    // Initialize the rasters.
    g_debug(" -- crowns: preparing rasters");
    Raster<float> inrast(config.crownsSmoothedCHM);
    Raster<uint32_t> outrast(config.crownsCrownsRaster, 1, inrast);
    outrast.nodata(0);
    outrast.fill(0);

    double nodata = inrast.nodata();

    if(m_callbacks)
        m_callbacks->stepCallback(0.02);
    
    // Initialize the database, get the treetop count and estimate the buffer size.
    g_debug(" -- crowns: preparing database")
    uint64_t geomCount;
    SQLite db(config.crownsTreetopsDatabase);
    db.getGeomCount(&geomCount);
    g_debug(" -- crowns: processing " << geomCount << " tree tops");

    if(m_callbacks)
        m_callbacks->stepCallback(0.03);
    
    // Build the list of offsets for D8 or D4 search.
    std::vector<std::pair<int, int> > offsets; // pairs of col, row
    offsets.push_back(std::make_pair(-1, -1));
    offsets.push_back(std::make_pair(-1, 0));
    offsets.push_back(std::make_pair(-1, 1));
    offsets.push_back(std::make_pair(0, -1));
    offsets.push_back(std::make_pair(0, 1));
    offsets.push_back(std::make_pair(1, -1));
    offsets.push_back(std::make_pair(1, 0));
    offsets.push_back(std::make_pair(1, 1));
    
    // The number of extra rows above and below the buffer.
    int32_t bufRows = (int) std::ceil(g_abs(config.crownsRadius / inrast.resolutionY()));
    // The height of the row, not including disposable buffer. Use bufRows as lower bound to
    // avoid read error later (keeps min row index to >=0)
    int32_t rowStep = g_max(bufRows, (int) g_abs(std::ceil((double) batchSize / geomCount * inrast.rows()) / inrast.resolutionY()));
    // The total height of the buffer
    int32_t rowHeight = rowStep + bufRows * 2;
    std::atomic<uint32_t> curRow(0);

    #pragma omp parallel
    {
        // To keep track of visited cells.
        std::vector<bool> visited((uint64_t) inrast.cols() * rowHeight);
        MemRaster<float> buf(inrast.cols(), rowHeight);
        MemRaster<uint32_t> blk(inrast.cols(), rowHeight);
        buf.fill(inrast.nodata());
        blk.fill(0);

        #pragma omp for
        for (int32_t row = 0; row < inrast.rows(); row += rowStep) {
            curRow += rowStep;

            // Load the tree tops for the strip.
            Bounds bounds(inrast.toX(0), inrast.toY(row - bufRows), inrast.toX(inrast.cols()), inrast.toY(row + rowStep + bufRows));
            std::vector<std::unique_ptr<Point> > tops;
            
            #pragma omp critical(__crowns_getpoints)
            db.getPoints(tops, bounds);

            #pragma omp critical(__crowns_readbuf) 
            inrast.readBlock(0, row == 0 ? row : row - bufRows, buf, 0, row == 0 ? bufRows : 0);

            // Convert the Tops to Nodes.
            std::queue<std::unique_ptr<Node> > q;
            for (const std::unique_ptr<Point> &t : tops) {
                int col = inrast.toCol(t->x);
                int row = inrast.toRow(t->y);
                int id = atoi(t->fields["id"].c_str());
                std::unique_ptr<Node> nd(new Node(id, col, row, t->z, col, row, t->z));
                q.push(std::move(nd));
            }

            // Run through the queue.
            while (q.size()) {
                std::unique_ptr<Node> n = std::move(q.front());
                q.pop();

                blk.set(n->c, n->r - row + bufRows, (uint32_t) n->id);

                for (const std::pair<int, int> &offset : offsets) {
                    int c = n->c + offset.first;
                    int r = n->r + offset.second;

                    if (r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols()) continue;
                    if (r - row + bufRows < 0 || r - row + bufRows >= buf.rows()) continue;

                    uint64_t idx = (uint64_t) (r - row + bufRows) * inrast.cols() + c;
                    if (visited[idx])
                        continue;

                    double v = buf.get(idx);
                    if (v != nodata                                                             // is not nodata
                            && v < n->z                                                         // is less than the neighbouring pixel
                            && v >= config.crownsMinHeight                                      // is greater than the min height
                            && (v / n->tz) >= config.crownsHeightFraction                       // is greater than the threshold height
                            && g_sq(n->tc - c) + g_sq(n->tr - r) <= g_sq(config.crownsRadius)   // is within the radius
                    ) {
                        std::unique_ptr<Node> nd(new Node(n->id, c, r, v, n->tc, n->tr, n->tz));
                        q.push(std::move(nd));
                        visited[idx] = true;
                    }
                }
            }

            if(m_callbacks)
                m_callbacks->stepCallback(0.03 + ((float) curRow / inrast.rows()) * 0.97);

            if (row > 0 && (row + bufRows) >= inrast.rows())
                continue;

            #pragma omp critical(__b)
            outrast.writeBlock(0, row, blk, 0, bufRows);
        }
    }
    
    if(m_callbacks)
        m_callbacks->stepCallback(1.0);

}

/*
 * treecrowns.cpp
 *
 *  Created on: May 3, 2016
 *      Author: rob
 */

#include <queue>
#include <iostream>
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

            /**
             * Represents a grid cell, and maintains some properties of the seed that originated it.
             */
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

            /**
             * Returns the distance between the given column/row pair in map units.
             */
            double dist(int tc, int tr, int c, int r, double resolution) {
                return std::sqrt(g_sq((double) (tc - c)) + g_sq((double) (tr - r))) * resolution;
            }

            /**
             * Returns true if the pixel at the center of the given
             * raster is the maximum value in the raster.
             */
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

            /**
             * Guess a value for a cell, based on its neighbours.
             */
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
x(x), y(y), z(z), uz(z),
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
    if(m_callbacks) {
        m_callbacks->fileCallback(0.25);
        m_callbacks->overallCallback(0.33);
    }
    Raster<float> in(config.smoothOriginalCHM);
    Raster<float> out(config.smoothSmoothedCHM, 1, in);
    in.smooth(out, config.smoothStdDev, config.smoothWindowSize);
    if(m_callbacks) {
        m_callbacks->fileCallback(1.0);
        m_callbacks->overallCallback(0.33);
    }
}

void Treetops::treetops(const TreetopsConfig &config) {
    if(m_callbacks) {
        m_callbacks->fileCallback(0.1);
        m_callbacks->overallCallback(0.66);
    }

    config.checkTops();

    // Initialize input rasters.
    g_debug(" -- tops: opening rasters");
    Raster<float> original(config.topsOriginalCHM);
    Raster<float> smoothed(config.topsSmoothedCHM);

    // Prepare database.
    g_debug(" -- tops: preparing db");
    std::map<std::string, int> fields;
    fields["id"] = 1;
    SQLite db(config.topsTreetopsDatabase, SQLite::POINT, config.srid, fields); // TODO: Get SRID from raster.
    db.makeFast();
    db.dropGeomIndex();
    db.setCacheSize(config.tableCacheSize);
    db.clear(); // TODO: Faster to delete it and start over.

    // This is the size of the cache used by each thread.
    uint32_t cachedRows = g_max(100, (config.rowCacheSize / original.rows() / sizeof (float))); // TODO: Make row cache configurable.
    uint32_t blockHeight = config.topsWindowSize * 2 + cachedRows;
    uint32_t numBlocks = original.rows() / blockHeight;
    uint64_t topId = 0;
    uint64_t topCount = 0;
    uint32_t curBlock = 0;
    
    #pragma omp parallel
    {
        g_debug(" -- tops: processing " << omp_get_thread_num());
        MemRaster<float> blk(original.cols(), blockHeight);
        blk.nodata(original.nodata());
        blk.fill(blk.nodata());

        std::map<uint64_t, std::unique_ptr<Top> > tops0;
        uint64_t topCount0 = 0;
        uint32_t curRow;
        double max;

        #pragma omp for
        for (uint32_t blockNum = 0; blockNum < numBlocks; ++blockNum) {
            curRow = blockNum * blockHeight - config.topsWindowSize;
            
            #pragma omp atomic
            curBlock++;
            
            #pragma omp critical(c)
            smoothed.readBlock(0, g_max(0, curRow), blk);

            for (uint32_t r = 0; r < blockHeight; ++r) {
                for (int c = 0; c < original.cols(); ++c) {

                    uint64_t id = ((uint64_t) c << 32) | (r + curRow);

                    if (blk.get(c, r) >= config.topsMinHeight &&
                            isMaxCenter(blk, c, r, config.topsWindowSize, &max)) {

                        // Get the original height from the unsmoothed raster.
                        // TODO: May not correspond to the smoothed top
                        double umax = original.get(c, r + curRow);

                        #pragma omp atomic
                        ++topId;

                        std::unique_ptr<Top> t(new Top(
                                topId,
                                original.toCentroidX(c + config.topsWindowSize / 2), // center of pixel
                                original.toCentroidY(r + curRow + config.topsWindowSize / 2),
                                max,
                                umax,
                                c + config.topsWindowSize / 2,
                                r + curRow + config.topsWindowSize / 2
                                ));

                        tops0[id] = std::move(t);
                        ++topCount0;
                    }
                }
            }
            
            if(m_callbacks) {
                m_callbacks->fileCallback((float) curBlock / numBlocks * 0.5);
                m_callbacks->overallCallback(0.66);
            }
            
        }

        #pragma omp atomic
        topCount += topCount0; // TODO: Doesn't work for status because there's no barrier here.

        uint64_t b = 0;
        uint64_t batch = db.maxAddPointCount();

        std::vector<std::unique_ptr<Point> > points;
        for (auto it = tops0.begin(); it != tops0.end(); ++it) {

            Top *t = it->second.get();

            std::map<std::string, std::string> fields;
            fields["id"] = std::to_string(t->id);

            Point *pt = new Point(t->x, t->y, t->uz, fields);
            points.push_back(std::unique_ptr<Point>(pt));

            if (++b % batch == 0 || b == topCount0) {
                #pragma omp critical(b)
                {
                    g_debug("inserting " << points.size() << " points");
                    db.begin();
                    db.addPoints(points);
                    db.commit();
                }
                points.clear();
                if(m_callbacks) {
                    m_callbacks->fileCallback((float) b / topCount * 0.5 + 0.5);
                    m_callbacks->overallCallback(0.66);
                }
            }
        }
    }

    if (config.buildIndex) {
        g_debug(" -- tops: building index");
        db.createGeomIndex();
    }
    db.makeSlow();

    g_debug(" -- tops: done.");
    if(m_callbacks) {
        m_callbacks->fileCallback(1.0);
        m_callbacks->overallCallback(0.66);
    }
}

void Treetops::treecrowns(const TreetopsConfig &config) {
    config.checkCrowns();

    if(m_callbacks) {
        m_callbacks->fileCallback(0.0);
        m_callbacks->overallCallback(0.99);
    }

    g_debug(" -- crowns: delineating crowns...");
    int batchSize = 1000000;

    // Initialize the rasters.
    g_debug(" -- crowns: preparing rasters");
    Raster<float> inrast(config.crownsSmoothedCHM);
    Raster<uint32_t> outrast(config.crownsCrownsRaster, 1, inrast);
    outrast.nodata(0);
    outrast.fill(0);

    double nodata = inrast.nodata();

    // Initialize the database, get the treetop count and estimate the buffer size.
    g_debug(" -- crowns: preparing database")
    SQLite db(config.crownsTreetopsDatabase);
    uint64_t geomCount;
    db.getGeomCount(&geomCount);
    g_debug(" -- crowns: processing " << geomCount << " tree tops");

    // The number of extra rows above and below the buffer.
    int bufRows = (int) std::ceil(g_abs(config.crownsRadius / inrast.resolutionY()));
    // The height of the row, not including disposable buffer. Use bufRows as lower bound to
    // avoid read error later (keeps min row index to >=0)
    int rowStep = g_max(bufRows, (int) g_abs(std::ceil((double) batchSize / geomCount * inrast.rows()) / inrast.resolutionY()));
    // The total height of the buffer
    int rowHeight = rowStep + bufRows * 2;
    int rowCompleted = 0;

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

    int curRow = 0;
    
    #pragma omp parallel for
    for (int row = 0; row < inrast.rows(); row += rowStep) {

        #pragma omp atomic
        curRow += rowStep;
        
        g_debug("Processing " << batchSize << " of " << geomCount << " points.");
        g_debug(" - row " << row << " of " << inrast.rows() << "; step: " << rowStep);

        // To keep track of visited cells.
        std::vector<bool> visited((uint64_t) inrast.cols() * rowHeight);
        MemRaster<float> buf(inrast.cols(), rowHeight);
        MemRaster<uint32_t> blk(inrast.cols(), rowHeight);
        buf.fill(inrast.nodata());
        blk.fill(0);

        // Load the tree tops for the strip.
        Bounds bounds(inrast.toX(0), inrast.toY(row - bufRows), inrast.toX(inrast.cols()), inrast.toY(row + rowStep + bufRows));
        std::vector<std::unique_ptr<Point> > tops;
        #pragma omp critical(crowns_getpoints)
        {
            db.getPoints(tops, bounds);
        }
        g_debug(" - reading from source: " << 0 << ", " << (row == 0 ? row : row - bufRows) << ", [buf], " << 0 << ", " << (row == 0 ? bufRows : 0));
        #pragma omp critical(crowns_readbuf) 
        {
            inrast.readBlock(0, row == 0 ? row : row - bufRows, buf, 0, row == 0 ? bufRows : 0);
        }

        // Convert the Tops to Nodes.
        std::queue<std::unique_ptr<Node> > q;
        for (const std::unique_ptr<Point> &t : tops) {
            int col = inrast.toCol(t->x);
            int row = inrast.toRow(t->y);
            int id = atoi(t->fields["id"].c_str());
            q.push(std::unique_ptr<Node>(new Node(id, col, row, t->z, col, row, t->z)));
        }

        // Run through the queue.
        while (q.size()) {

            std::unique_ptr<Node> n = std::move(q.front());
            q.pop();

            blk.set(n->c, n->r - row + bufRows, (uint32_t) n->id);

            for (const std::pair<int, int> &offset : offsets) {
                int c = n->c + offset.first;
                int r = n->r + offset.second;

                //g_debug(" -- " << n->c << ", " << n->r << "; " << c << ", " << r);
                if (r < 0 || c < 0 || r >= inrast.rows() || c >= inrast.cols()) continue;
                if (r - row + bufRows < 0 || r - row + bufRows >= buf.rows()) continue;

                uint64_t idx = (uint64_t) (r - row + bufRows) * inrast.cols() + c;
                if (visited[idx])
                    continue;

                double v = buf.get(idx);
                if (v != nodata // is not nodata
                        && v < n->z // is less than the neighbouring pixel
                        && v >= config.crownsMinHeight // is greater than the min height
                        && (v / n->tz) >= config.crownsHeightFraction // is greater than the threshold height
                        && std::pow(n->tc - c, 2) + std::pow(n->tr - r, 2) <= std::pow(config.crownsRadius, 2) // is within the radius
                ) {
                    q.push(std::unique_ptr<Node>(new Node(n->id, c, r, v, n->tc, n->tr, n->tz)));
                    visited[idx] = true;
                }
            }
        }

        if(m_callbacks) {
            m_callbacks->fileCallback((float) curRow / inrast.rows());
            m_callbacks->overallCallback(0.66);
        }

        g_debug(" - tmp block: cols: " << blk.cols() << ", rows: " << (row == 0 ? rowStep + bufRows : rowStep));
        if (row > 0 && (row + bufRows) >= inrast.rows())
            continue;
        #pragma omp critical(b)
        {
            outrast.writeBlock(0, row, blk, 0, bufRows);
            rowCompleted += rowStep;
        }
    }

    if(m_callbacks) {
        m_callbacks->fileCallback(1.0);
        m_callbacks->overallCallback(1.0);
    }
}

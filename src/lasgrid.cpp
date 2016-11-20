/*
 * Grids a point cloud represented by one or more LAS files.
 * Can produce grids of from intensity and elevation, using
 * minimum, maximum, mean, std dev, density, variance and count.
 *
 * Authored by: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <cstring>
#include <math.h>
#include <exception>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "lasgrid.hpp"
#include "lasutil.hpp"

namespace fs = boost::filesystem;
namespace alg = boost::algorithm;

using namespace geotools::util;
using namespace geotools::raster;

namespace geotools {

    namespace las {

        namespace lasgrid_config {

            double defaultResolution = 2.0;
            double defaultRadius = std::sqrt(g_sq(defaultResolution / 2.0) * 2.0);
            bool defaultSnapToGrid = true;
            int defaultType = TYPE_MEAN;
            int defaultAttribute = ATT_HEIGHT;
            unsigned char defaultAngleLimit = 180;
            std::set<int> defaultClasses = {2};
            std::map<std::string, int> types = {
                {"Minimum", TYPE_MIN},
                {"Maximum", TYPE_MAX},
                {"Mean", TYPE_MEAN},
                {"Density", TYPE_DENSITY},
                {"Sample Variance", TYPE_VARIANCE},
                {"Sample Std. Dev.", TYPE_STDDEV},
                {"Population Variance", TYPE_PVARIANCE},
                {"Population Std. Dev.", TYPE_PSTDDEV},
                {"Count", TYPE_COUNT},
                {"Quantile", TYPE_QUANTILE},
                {"Median", TYPE_MEDIAN}
            };
            std::map<std::string, int> attributes = {
                {"Height", ATT_HEIGHT},
                {"Intensity", ATT_INTENSITY}
            };

        } // config

        namespace lasgrid_util {

            /**
             * Interpret the value of a string attribute name, return the constant int value.
             */
            int parseAtt(char *attStr) {
                if (!strcmp("intensity", attStr)) {
                    return ATT_INTENSITY;
                } else if (!strcmp("height", attStr)) {
                    return ATT_HEIGHT;
                }
                return 0;
            }

            /**
             * Interpret the output type and return the constant int value.
             */
            int parseType(char *typeStr) {
                if (!strcmp("min", typeStr)) {
                    return TYPE_MIN;
                } else if (!strcmp("max", typeStr)) {
                    return TYPE_MAX;
                } else if (!strcmp("mean", typeStr)) {
                    return TYPE_MEAN;
                } else if (!strcmp("density", typeStr)) {
                    return TYPE_DENSITY;
                } else if (!strcmp("variance", typeStr)) {
                    return TYPE_VARIANCE;
                } else if (!strcmp("stddev", typeStr)) {
                    return TYPE_STDDEV;
                } else if (!strcmp("pvariance", typeStr)) {
                    return TYPE_PVARIANCE;
                } else if (!strcmp("pstddev", typeStr)) {
                    return TYPE_PSTDDEV;
                } else if (!strcmp("count", typeStr)) {
                    return TYPE_COUNT;
                } else if (!strcmp("median", typeStr)) {
                    return TYPE_MEDIAN;
                }
                return 0;
            }

            /**
             * Comparator for sorting doubles.
             */
            int _fcmp(const void * a, const void * b) {
                const double * aa = (const double *) a;
                const double * bb = (const double *) b;
                return (*aa > *bb) - (*bb > *aa);
            }

            void vector_dealloc(std::vector<double> *item) {
                delete item;
            }

            /**
             * Returns true if the point is within the radius associated with a cell's centroid.
             * @param px The x coordinate of the point.
             * @param py The y coordinate of the point.
             * @param col The column of the cell of interest.
             * @param row The row of the cell of interest.
             * @param radius The radius around the cell's centroid.
             * @param resolution The resolution of the output raster.
             * @param bounds The bounds of the raster.
             */
            bool inRadius(double px, double py, int col, int row, double radius,
                    double resolution, Bounds &bounds) {
                if (radius == 0.0) return true;
                // If a radius is given, extract the x and y of the current cell's centroid
                // and measure its distance (squared) from the point.
                double x = col * resolution + bounds[0] + resolution * 0.5;
                double y = row * resolution + bounds[1] + resolution * 0.5;
                // If the cell is outside the radius, ignore it.
                double r = sqrt(g_sq(x - px) + g_sq(y - py));
                return r <= radius;
            }

        } // util

        void LasGrid::setCallbacks(geotools::util::Callbacks *callbacks) {
            m_callbacks = callbacks;
        }

        void LasGrid::lasgrid(std::string &dstFile, std::vector<std::string> &files, std::set<int> &classes,
                int crs, int attribute, int type, double radius,
                double resolution, Bounds &bounds, unsigned char angleLimit, bool fill, bool snap) {

            if (resolution <= 0.0)
                g_argerr("Resolution must be > 0.");
            if (radius <= 0.0) {
                radius = std::sqrt(g_sq(resolution / 2.0) * 2.0);
                g_warn("Radius invalid; using " << radius);
            }
            if (files.size() == 0)
                g_argerr("At least one input file is required.");
            if (dstFile.empty())
                g_argerr("An output file is required.");
            if (attribute <= 0)
                g_argerr("An attribute is required.");
            if (type <= 0)
                g_argerr("A valid type is required.");
            if (classes.size() == 0)
                g_warn("No classes given. Matching all classes.");
            if (angleLimit <= 0)
                g_argerr("Angle limit must be greater than zero.");

            g_debug("Radius: " << radius);
            g_debug("Resolution: " << resolution);
            g_debug("Files: " << files.size());
            g_debug("Destination: " << dstFile);
            g_debug("Attribute: " << attribute);
            g_debug("Type: " << type);
            g_debug("Classes: " << classes.size());
            g_debug("Angle Limit: " << angleLimit);

            using namespace geotools::las::lasgrid_util;

            MemRaster<double> grid1;
            MemRaster<int> counts;
            MemRaster<std::vector<double>* > qGrid;

            liblas::ReaderFactory rf;
            std::vector<unsigned int> indices;
            Bounds bounds1;
            bounds1.collapse();
            g_debug("Total bounds initial: " << bounds1.print());

            for (unsigned int i = 0; i < files.size(); ++i) {
                g_debug("Checking file " << files[i]);
                std::ifstream in(files[i].c_str(), std::ios::in | std::ios::binary);
                liblas::Reader r = rf.CreateWithStream(in);
                liblas::Header h = r.GetHeader();
                Bounds bounds0;
                if (!LasUtil::computeLasBounds(h, bounds0, 2))
                    LasUtil::computeLasBounds(r, bounds0, 2); // If the header bounds are bogus.
                g_debug("File bounds " << files[i] << ": " << bounds0.print());
                in.close();
                if (bounds.intersects(bounds0, 2)) {
                    indices.push_back(i);
                    bounds1.extend(bounds0);
                    g_debug("Total bounds: " << bounds1.print());
                }
            }

            bounds.collapse();
            bounds.extend(bounds1);

            if (snap)
                bounds.snap(resolution);

            g_debug(bounds.print());

            // Prepare grid
            int cols = bounds.cols(resolution);
            int rows = bounds.rows(resolution);

            g_debug("Raster size: " << cols << ", " << rows << "; Cell radius: " << radius);

            // For types other than count, we need a double grid to manage sums.
            if (type != TYPE_COUNT) {
                grid1.init(cols, rows);
                grid1.fill(-9999.0);
                g_debug("Created grid1: " << cols << ", " << rows);
            }

            // For the median grid.
            switch (type) {
                case TYPE_VARIANCE:
                case TYPE_STDDEV:
                case TYPE_PVARIANCE:
                case TYPE_PSTDDEV:
                case TYPE_QUANTILE:
                case TYPE_MEDIAN:
                    qGrid.init(cols, rows);
                    qGrid.setDeallocator(&vector_dealloc);
                    for (size_t i = 0; i < qGrid.size(); ++i)
                        qGrid.set(i, new std::vector<double>());
                    g_debug("Created qGrid: " << cols << ", " << rows);
                    break;
            }

            // Create a grid for maintaining counts.
            g_debug("Creating counts grid: " << cols << ", " << rows);
            counts.init(cols, rows);
            counts.fill(0);

            g_debug("Using " << indices.size() << " of " << files.size() << " files.");

            // Process files
            for (unsigned int j = 0; j < indices.size(); ++j) {

                unsigned int i = indices[j];
                std::ifstream in(files[i].c_str());
                liblas::Reader reader = rf.CreateWithStream(in);
                liblas::Header header = reader.GetHeader();

                if (m_callbacks)
                    m_callbacks->overallCallback(((j + 0.5f) / indices.size()) * 0.75f);

                size_t curPt = 0;
                size_t numPts = header.GetPointRecordsCount();

                while (reader.ReadNextPoint()) {
                    liblas::Point pt = reader.GetPoint();

                    if (curPt % 1000 == 0) {
                        if (m_callbacks)
                            m_callbacks->fileCallback((curPt + 1.0) / numPts);
                    }
                    ++curPt;

                    if (g_abs(pt.GetScanAngleRank()) > angleLimit)
                        continue;

                    double px = pt.GetX();
                    double py = pt.GetY();

                    // Check if in bounds, but only if clipping is desired.
                    if (!bounds.contains(px, py))
                        continue;
                    // If this point is not in the class list, skip it.
                    if (!Util::inList(classes, pt.GetClassification().GetClass()))
                        continue;

                    // Get either the height or intensity value.
                    double pz;
                    if (attribute == ATT_INTENSITY) {
                        pz = pt.GetIntensity();
                    } else { // ATT_HEIGHT
                        pz = pt.GetZ();
                    }

                    // Convert x and y, to col and row.
                    int c = (int) ((px - bounds.minx()) / resolution);
                    int r = (int) ((py - bounds.miny()) / resolution);

                    // If the radius is > 0, compute the size of the window.
                    int offset = (int) (radius * 2) / resolution;
                    for (int cc = g_max(0, c - offset); cc < g_min(cols, c + offset + 1); ++cc) {
                        for (int rr = g_max(0, r - offset); rr < g_min(rows, r + offset + 1); ++rr) {
                            // If the coordinate is out of the cell's radius, continue.
                            if (!inRadius(px, py, cc, rr, radius, resolution, bounds))
                                continue;
                            // Compute the grid index. The rows are assigned from the bottom.
                            int idx = (rows - rr - 1) * cols + cc;
                            //g_debug("idx: " << idx);
                            counts.set(idx, counts.get(idx) + 1);

                            switch (type) {
                                case TYPE_MIN:
                                    if (counts[idx] == 1 || pz < grid1[idx])
                                        grid1.set(idx, pz);
                                    break;
                                case TYPE_MAX:
                                    if (counts[idx] == 1 || pz > grid1[idx])
                                        grid1.set(idx, pz);
                                    break;
                                case TYPE_MEAN:
                                    if (counts[idx] == 1) {
                                        grid1.set(idx, pz);
                                    } else {
                                        grid1.set(idx, grid1.get(idx) + pz);
                                    }
                                    break;
                                case TYPE_VARIANCE:
                                case TYPE_STDDEV:
                                case TYPE_PVARIANCE:
                                case TYPE_PSTDDEV:
                                case TYPE_QUANTILE:
                                case TYPE_MEDIAN:
                                    qGrid[idx]->push_back(pz);
                                    break;
                            }
                        }
                    }
                }

                if (m_callbacks) {
                    m_callbacks->fileCallback(1.0f);
                    m_callbacks->overallCallback(((j + 1.0f) / indices.size()) * 0.75f);
                }

            }

            if (m_callbacks)
                m_callbacks->overallCallback(0.75f);

            // Calculate cells or set nodata.
            // Welford's method for variance.
            switch (type) {
                case TYPE_MEAN:
                    for (size_t i = 0; i < (size_t) cols * rows; ++i) {
                        if (counts[i] > 0)
                            grid1.set(i, grid1.get(i) / counts.get(i));
                    }
                    break;
                case TYPE_PVARIANCE:
                    for (size_t i = 0; i < (size_t) cols * rows; ++i) {
                        if (counts[i] > 1) {
                            double m = 0;
                            double s = 0;
                            int k = 1;
                            for (unsigned int j = 0; j < qGrid[i]->size(); ++j) {
                                double v = qGrid[i]->at(j);
                                double oldm = m;
                                m = m + (v - m) / k;
                                s = s + (v - m) * (v - oldm);
                                ++k;
                            }
                            grid1.set(i, s / qGrid[i]->size());
                        } else {
                            grid1.set(i, 0);
                        }
                    }
                    break;
                case TYPE_VARIANCE:
                    for (size_t i = 0; i < (size_t) cols * rows; ++i) {
                        if (counts[i] > 1) {
                            double m = 0;
                            double s = 0;
                            int k = 1;
                            for (unsigned int j = 0; j < qGrid[i]->size(); ++j) {
                                double v = qGrid[i]->at(j);
                                double oldm = m;
                                m = m + (v - m) / k;
                                s = s + (v - m) * (v - oldm);
                                ++k;
                            }
                            grid1.set(i, s / (qGrid[i]->size() - 1));
                        } else {
                            grid1.set(i, 0);
                        }
                    }
                    break;
                case TYPE_PSTDDEV:
                    for (size_t i = 0; i < (size_t) cols * rows; ++i) {
                        if (counts[i] > 1) {
                            double m = 0;
                            double s = 0;
                            int k = 1;
                            for (unsigned int j = 0; j < qGrid[i]->size(); ++j) {
                                double v = qGrid[i]->at(j);
                                double oldm = m;
                                m = m + (v - m) / k;
                                s = s + (v - m) * (v - oldm);
                                ++k;
                            }
                            grid1.set(i, std::sqrt(s / qGrid[i]->size()));
                        } else {
                            grid1.set(i, 0);
                        }
                    }
                    break;
                case TYPE_STDDEV:
                    for (size_t i = 0; i < (size_t) cols * rows; ++i) {
                        if (counts[i] > 1) {
                            double m = 0;
                            double s = 0;
                            int k = 1;
                            for (unsigned int j = 0; j < qGrid[i]->size(); ++j) {
                                double v = qGrid[i]->at(j);
                                double oldm = m;
                                m = m + (v - m) / k;
                                s = s + (v - m) * (v - oldm);
                                ++k;
                            }
                            grid1.set(i, std::sqrt(s / (qGrid[i]->size() - 1)));
                        } else {
                            grid1.set(i, 0);
                        }
                    }
                    break;
                case TYPE_DENSITY:
                {
                    double r2 = g_sq(resolution);
                    for (size_t i = 0; i < (size_t) cols * rows; ++i) {
                        if (counts[i] > 0) {
                            grid1.set(i, (double) counts[i] / r2);
                        } else {
                            grid1.set(i, 0.0);
                        }
                    }
                }
                    break;
                case TYPE_MEDIAN:
                    for (size_t i = 0; i < (size_t) cols * rows; ++i) {
                        if (counts[i] > 0) {
                            std::sort(qGrid[i]->begin(), qGrid[i]->end());
                            int size = qGrid[i]->size();
                            if (size % 2 == 0) {
                                int idx = size / 2;
                                grid1.set(i, (qGrid[i]->at(idx - 1) + qGrid[i]->at(idx)) / 2.0);
                            } else {
                                grid1.set(i, (qGrid[i]->at(size / 2.0)));
                            }
                        }
                    }
                    break;
            }

            if (m_callbacks)
                m_callbacks->overallCallback(0.85f);

            if (type == TYPE_COUNT) {
                // TODO: Determine resolution sign properly.
                Raster<int> rast(dstFile, 1, bounds, resolution, -resolution, -1, crs);
                rast.writeBlock(counts);
            } else {
                Raster<double> rast(dstFile, 1, bounds, resolution, -resolution, -9999.0, crs);
                rast.writeBlock(grid1);
                //if(fill)
                //	rast.voidFillIDW(resolution);
            }

            if (m_callbacks)
                m_callbacks->overallCallback(1.0f);

        }

    } // las

} // geotools


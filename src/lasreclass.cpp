/*
 * This tool reassigns the classes on each point.
 *
 *  Author: Rob Skelly rob@dijital.ca
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <list>
#include <climits>
#include <memory>
#include <iomanip>
#include <queue>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <liblas/liblas.hpp>

#include "geotools.h"

#define TIME_GAP 50.0

namespace geotools {

    namespace las {

        namespace util {

            static int __seg_id = 0;

            class Seg {
            public:
                int id = ++__seg_id;
                double start = nan("");
                double end = nan("");
                liblas::Header *header = nullptr;
                liblas::Writer *writer = nullptr;
                std::ofstream out;
                bool inited = false;
                double bounds[6] = {G_DBL_MAX_POS, G_DBL_MAX_POS, G_DBL_MAX_NEG, G_DBL_MAX_NEG, G_DBL_MAX_POS, G_DBL_MAX_NEG};
                bool marked;
                int ptCount = 0;
                int ptRecByRetCount[5] = {0, 0, 0, 0, 0};

                Seg(double start, double end) {
                    this->start = start;
                    this->end = end;
                }

                void initialize(const std::string &filename, liblas::Header &h) {
                    std::string file = filename + "_" + std::to_string(id) + ".las";
                    g_trace("Seg initialized as " << file);
                    out.open(file.c_str(), std::ios::out | std::ios::binary);
                    header = new liblas::Header(h);
                    writer = new liblas::Writer(out, *header);
                    inited = true;
                }

                ~Seg() {
                    finalize();
                }

                void addPoint(liblas::Point &pt) {
                    if (pt.GetX() < bounds[0]) bounds[0] = pt.GetX();
                    if (pt.GetY() < bounds[1]) bounds[1] = pt.GetY();
                    if (pt.GetX() > bounds[2]) bounds[2] = pt.GetX();
                    if (pt.GetY() > bounds[3]) bounds[3] = pt.GetY();
                    if (pt.GetZ() < bounds[4]) bounds[4] = pt.GetZ();
                    if (pt.GetZ() > bounds[5]) bounds[5] = pt.GetZ();
                    pt.SetPointSourceID(id);
                    writer->WritePoint(pt);
                    ++ptCount;
                    ++ptRecByRetCount[pt.GetReturnNumber() - 1];
                }

                void finalize() {
                    if (header != nullptr) {
                        header->SetMin(bounds[0], bounds[1], bounds[4]);
                        header->SetMax(bounds[2], bounds[3], bounds[5]);
                        header->SetPointRecordsCount(ptCount);
                        for (size_t i = 0; i < 5; ++i)
                            header->SetPointRecordsByReturnCount(i + 1, ptRecByRetCount[i]);
                        writer->SetHeader(*header);
                        writer->WriteHeader();
                        out.close();
                        delete writer;
                        delete header;
                        header = nullptr;
                        writer = nullptr;
                    }
                }

                bool intersects(Seg *seg) {
                    return !(seg->end < start || seg->start > end);
                }

                bool near(Seg *seg) {
                    return ((start - seg->end) < 1.0 && start > seg->end) || ((seg->start - end) < 1.0 && seg->start > end);
                }

                bool insert(Seg *seg) {
                    g_trace("Insert: " << seg->id << ": " << seg->start << " > " << seg->end);
                    if (seg->id == id) {
                        g_trace("> same seg.");
                        return false;
                    } else if (intersects(seg) || near(seg)) {
                        g_trace("> " << id << " joins " << seg->id);
                        start = g_min(start, seg->start);
                        end = g_max(end, seg->end);
                    } else {
                        g_trace("> " << id << " no relationship " << seg->id);
                        return false;
                    }
                    g_trace("New span: " << id << ": " << start << " > ");
                    return true;
                }
            };

            struct {

                bool operator()(Seg *a, Seg *b) {
                    return a->start < b->start;
                }
            } segsort;

            void normalizeFlightLines(std::vector<Seg*> &flightLines) {
                int count = flightLines.size();
                std::sort(flightLines.begin(), flightLines.end(), segsort);
                std::set<Seg*> output;
                Seg *seg = flightLines[0];
                output.insert(seg);
                for (size_t i = 1; i < flightLines.size(); ++i) {
                    if (!seg->insert(flightLines[i])) {
                        seg = flightLines[i];
                        output.insert(seg);
                    } else {
                        delete flightLines[i];
                    }
                }
                flightLines.resize(0);
                flightLines.assign(output.begin(), output.end());
                g_trace(flightLines.size() << " remaining of " << count);
                std::sort(flightLines.begin(), flightLines.end(), segsort);
                for (Seg *seg : flightLines)
                    g_trace(" -- " << seg->id << ": " << seg->start << " > " << seg->end);
            }

            // TODO: Replace with interval tree.

            int findFlightLine(std::vector<Seg*> &flightLines, double time) {
                for (Seg *seg : flightLines) {
                    if (time >= seg->start && time <= seg->end)
                        return seg->id;
                }
                g_trace("Seg for time not found: " << time);
                return 0;
            }

            void computeDirection(std::list<liblas::Point> &q, double *pdir) {
                liblas::Point pt0 = q.front();
                liblas::Point pt1 = q.back();
                double d = atan2(pt1.GetY() - pt0.GetY(), pt1.GetX() - pt0.GetX());
                while (d < 0)
                    d += G_PI * 2.0;
                while (d > G_PI * 2.0)
                    d -= G_PI * 2.0;
                *pdir = d;
            }

            void dumpQueue(std::list<liblas::Point> &pq, liblas::Writer &w, int cls = -1) {
                for (liblas::Point &p : pq) {
                    if (cls > -1) {
                        liblas::Classification c = p.GetClassification();
                        c.SetClass(cls);
                        p.SetClassification(c);
                    }
                    w.WritePoint(p);
                }
                pq.clear();
            }

            void recoverEdges(std::vector<std::string> &files, std::string &outfile) {

                if (outfile.size() == 0)
                    g_argerr("An output directory is required.");
                if (files.size() == 0)
                    g_argerr("At least one input file is required.");

                /* Loop over files and figure out which ones are relevant. */
                liblas::ReaderFactory rf;

                for (unsigned int i = 0; i < files.size(); ++i) {

                    std::string path = files[i];
                    std::string base = path.substr(path.find_last_of("/") + 1);
                    std::string newpath = outfile + "/" + base;

                    g_trace("Saving " << path << " to " << newpath);

                    std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
                    liblas::Reader r = rf.CreateWithStream(in);
                    liblas::Header h = r.GetHeader();

                    std::ofstream out(newpath.c_str(), std::ios::out | std::ios::binary);
                    liblas::Header wh(h);
                    liblas::Writer w(out, wh);

                    double dir0 = 0, dir1 = 0; // queue direction
                    std::list<liblas::Point> pq0;
                    std::list<liblas::Point> pq1;
                    int limit = 20;
                    double lastTime = 0.0;

                    while (r.ReadNextPoint()) {
                        //g_trace("Q0: " << pq0.size() << "; Q1: " << pq1.size());
                        liblas::Point pt = r.GetPoint();
                        //_print("time " << pt.GetTime());
                        double t = pt.GetTime();
                        if (lastTime == 0.0) {
                            lastTime = t;
                        } else if (t - lastTime > 0.4) {
                            lastTime = t;
                            dumpQueue(pq0, w);
                            dumpQueue(pq1, w);
                        }
                        pq0.push_back(pt);
                        if (pq0.size() <= limit) {
                            continue;
                        } else {
                            pq1.push_back(pq0.front());
                            pq0.pop_front();
                        }
                        if (pq1.size() == limit) {
                            computeDirection(pq0, &dir0);
                            computeDirection(pq1, &dir1);
                            //g_trace("Dirs: " << dir0 << ", " << dir1);
                            if (std::abs(std::abs(dir0) - std::abs(dir1)) > G_PI * 0.75) {
                                dumpQueue(pq0, w, 31);
                                dumpQueue(pq1, w, 31);
                                //g_trace("Flip " << dir0 << "," << dir1 << ": " << (int) pq1.front()->GetClassification().GetClass());
                            } else {
                                w.WritePoint(pq1.front());
                                pq1.pop_front();
                            }
                        }
                    }

                    dumpQueue(pq0, w);
                    dumpQueue(pq1, w);

                    in.close();
                    out.close();

                }
            }


        } // Util

        void mapClasses(std::vector<std::string> &files, std::string &outfile, std::map<int, int> &mappings) {

            if (mappings.size() == 0)
                g_argerr("At least one mapping is required.");
            g_trace("Mappings:");
            for (std::map<int, int>::iterator it = mappings.begin(); it != mappings.end(); ++it)
                g_trace(" " << it->first << " > " << it->second);
            if (outfile.size() == 0)
                g_argerr("An output directory (-o) is required.");
            if (files.size() == 0)
                g_argerr("At least one input file is required.");

            using namespace geotools::las::util;

            /* Loop over files and figure out which ones are relevant. */
            liblas::ReaderFactory rf;

            for (unsigned int i = 0; i < files.size(); ++i) {

                std::string path = files[i];
                std::string base = path.substr(path.find_last_of("/") + 1);
                std::string newpath = outfile + "/" + base;

                g_trace("Saving " << path << " to " << newpath);

                std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
                liblas::Reader r = rf.CreateWithStream(in);
                liblas::Header h = r.GetHeader();

                std::ofstream out(newpath.c_str(), std::ios::out | std::ios::binary);
                liblas::Header wh(h);
                liblas::Writer w(out, wh);

                while (r.ReadNextPoint()) {
                    liblas::Point pt = r.GetPoint();
                    liblas::Classification cl = pt.GetClassification();
                    if (mappings.find(cl.GetClass()) != mappings.end()) {
                        cl.SetClass(mappings[cl.GetClass()]);
                        pt.SetClassification(cl);
                    }
                    w.WritePoint(pt);
                }

                in.close();
                out.close();

            }
        }

        void recoverFlightlines(std::vector<std::string> &files, std::string &outfile, double timeGap, bool exportLas) {

            if (outfile.size() == 0)
                g_argerr("An output directory is required.");
            if (files.size() == 0)
                g_argerr("At least one input file is required.");
            if (timeGap <= 0)
                g_argerr("Time gap must be larger than zero.");

            using namespace geotools::las::util;

            /* Loop over files and figure out which ones are relevant. */
            liblas::ReaderFactory rf;

            // Store the IDs with the start or end time for that ID.
            // Subsequent files will search the table to find out if they have a 
            // continuation of this line.
            std::vector<Seg*> flightLines;
            double startTime = 0, endTime = 0;

            for (unsigned int i = 0; i < files.size(); ++i) {
                std::string path = files[i];

                g_trace("Checking " << path);

                std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
                liblas::Reader r = rf.CreateWithStream(in);
                liblas::Header h = r.GetHeader();

                startTime = 0;

                while (r.ReadNextPoint()) {
                    liblas::Point pt = r.GetPoint();
                    double time = pt.GetTime();
                    if (startTime == 0) {
                        startTime = endTime = time;
                    } else {
                        double gap = time - endTime;
                        if (gap < 0.0 || gap > timeGap) { // One what?
                            flightLines.push_back(new Seg(startTime, endTime));
                            g_trace("Time gap. Flightline: " << gap << "; " << startTime << " > " << endTime);
                            startTime = time;
                        }
                        endTime = time;
                    }
                }

                if (endTime != startTime) { // TODO: What is the time scale?)
                    flightLines.push_back(new Seg(startTime, endTime));
                    g_trace("Time gap. Flightline: " << startTime << " > " << endTime);
                }
                in.close();
            }

            normalizeFlightLines(flightLines);


            bool initSegs = true;
            std::map<int, Seg*> flMap;

            for (unsigned int i = 0; i < files.size(); ++i) {

                std::string path = files[i];

                g_trace("Processing " << path);

                std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
                liblas::Reader r = rf.CreateWithStream(in);
                liblas::Header h = r.GetHeader();

                if (initSegs) {
                    std::string base(outfile + "/flight_line");
                    for (Seg *seg : flightLines) {
                        flMap[seg->id] = seg;
                        seg->initialize(base, h);
                    }
                    initSegs = false;
                }

                while (r.ReadNextPoint()) {
                    liblas::Point pt = r.GetPoint();
                    Seg *seg = flMap[findFlightLine(flightLines, pt.GetTime())];
                    seg->addPoint(pt);
                }

                in.close();

            }

            for (Seg *seg : flightLines) {
                seg->finalize();
                delete seg;
            }

        }

        void recoverEdgesWEdgeMarker(std::vector<std::string> &files, std::string &outfile) {

            if (outfile.size() == 0)
                g_argerr("An output directory is required.");
            if (files.size() == 0)
                g_argerr("At least one input file is required.");

            using namespace geotools::las::util;

            /* Loop over files and figure out which ones are relevant. */
            liblas::ReaderFactory rf;

            for (unsigned int i = 0; i < files.size(); ++i) {

                std::string path = files[i];
                std::string base = path.substr(path.find_last_of("/") + 1);
                std::string newpath = outfile + "/" + base;

                g_trace("Saving " << path << " to " << newpath);

                std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
                liblas::Reader r = rf.CreateWithStream(in);
                liblas::Header h = r.GetHeader();

                std::ofstream out(newpath.c_str(), std::ios::out | std::ios::binary);
                liblas::Header wh(h);
                liblas::Writer w(out, wh);

                std::list<liblas::Point> pq0;
                std::list<liblas::Point> pq1;
                int limit = 20;

                while (r.ReadNextPoint()) {
                    //g_trace("Q0: " << pq0.size() << "; Q1: " << pq1.size());
                    liblas::Point pt = r.GetPoint();
                    //_print("time " << pt.GetTime());
                    pq0.push_back(pt);
                    if (pq0.size() > limit) {
                        // If the first queue is over the limit:
                        liblas::Point pt0 = pq0.front();
                        if (pt0.GetFlightLineEdge()) {
                            // If it is an edge, dump all the points on either side.
                            dumpQueue(pq0, w, 18);
                            dumpQueue(pq1, w, 18);
                        } else {
                            // Otherwise add it to the second queue
                            pq0.pop_front();
                            pq1.push_back(pt0);
                        }
                    }
                    if (pq1.size() == limit) {
                        // If the second queue is at the limit, write the last point.
                        w.WritePoint(pq1.front());
                        pq1.pop_front();
                    }
                }

                // TODO: If there's an edge in here, deal with it.
                dumpQueue(pq0, w);
                dumpQueue(pq1, w);

                in.close();
                out.close();

            }
        }

    } // las

} // geotools

void usage() {
    std::cerr << "Usage: lasreclass [options] srcfiles*\n"
            << " Performs certain operations on LAS attributes.\n"
            << " -o              Speficy an output folder (must exist).\n"
            << " -m <from> <to>  Specify a class mapping. Points from the first class are reclassified to the second class. Repeatable.\n"
            << " -f              Attempt to reconstruct flight lines using point times.\n"
            << " -fe             Export flightlines as single LAS files. Will have original filename with ID appended.\n"
            << " -fg             The time gap required to signify a break in flight lines. Defaul " << TIME_GAP << "\n";
}

int main(int argc, char ** argv) {

    bool mapping = false;
    bool flightlines = false;
    bool flExport = false;
    bool edges = false;
    double flTimeGap = TIME_GAP;
    std::string outfile;
    std::vector<std::string> files;
    std::map<int, int> mappings;

    /* Parse and check input. */
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "-m") {
            int from = atoi(argv[++i]);
            int to = atoi(argv[++i]);
            mappings[from] = to;
            mapping = true;
        } else if (arg == "-f") {
            flightlines = true;
        } else if (arg == "-fg") {
            flTimeGap = atof(argv[++i]);
        } else if (arg == "-fe") {
            flExport = true;
        } else if (arg == "-o") {
            outfile += argv[++i];
        } else if (arg == "-e") {
            edges = true;
        } else if (arg == "-v") {
            g_loglevel(1);
        } else {
            files.push_back(std::string(argv[i]));
        }
    }

    try {

        if (mapping) {

            geotools::las::mapClasses(files, outfile, mappings);

        } else if (flightlines) {

            geotools::las::recoverFlightlines(files, outfile, flTimeGap, flExport);

        } else if (edges) {

            // Only works with edge markers. 
            // TODO: Add other mechanisms of finding edges.
            geotools::las::recoverEdgesWEdgeMarker(files, outfile);

        } else {
            g_argerr("No command given.");
        }

    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        usage();
        return 1;
    }

    return 0;

}


#include <boost/filesystem.hpp>

#include <CGAL/Plane_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <liblas/liblas.hpp>

#include "pointnormalize.hpp"
#include "raster.hpp"
#include "util.hpp"
#include "lasreader.hpp"

using namespace geotools::point;
using namespace geotools::raster;
using namespace geotools::las;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::All_faces_iterator All_faces_iterator;
typedef Delaunay::Face Face;
typedef Delaunay::Face_handle Face_handle;

double computeArea(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
    double side0 = std::sqrt(std::pow(x1 - x2, 2.0) + std::pow(y1 - y2, 2.0) + std::pow(z1 - z2, 2.0));
    double side1 = std::sqrt(std::pow(x2 - x3, 2.0) + std::pow(y2 - y3, 2.0) + std::pow(z2 - z3, 2.0));
    double side2 = std::sqrt(std::pow(x3 - x1, 2.0) + std::pow(y3 - y1, 2.0) + std::pow(z3 - z1, 2.0));
    double s = (side0 + side1 + side2) / 2.0;
    return std::sqrt(s * (s - side0) * (s - side1) * (s - side2));
}

class FileSorter {
            private:
                double m_colSize;
                double m_rowSize;
            public:

                FileSorter(double colSize, double rowSize) :
                m_colSize(colSize), m_rowSize(rowSize) {
                }

                bool operator()(const std::string &a, const std::string &b) {
                    LASReader ar(a);
                    LASReader br(b);
                    Bounds ab = ar.bounds();
                    Bounds bb = br.bounds();
                    int idxa = ((int) (ab.miny() / m_rowSize)) * ((int) (ab.width() / m_colSize)) + ((int) (ab.minx() / m_colSize));
                    int idxb = ((int) (bb.miny() / m_rowSize)) * ((int) (bb.width() / m_colSize)) + ((int) (bb.minx() / m_colSize));
                    return idxa < idxb;
                }
            };
            
void PointNormalize::normalize(const PointNormalizeConfig &config, const Callbacks *callbacks) {

    if (config.pointOutputDir.empty())
        g_argerr("A point output dir must be provided.");
    
    if (!Util::mkdir(config.pointOutputDir))
        g_argerr("Couldn't create output dir " << config.pointOutputDir);

    using namespace boost::filesystem;

    path outDir(config.pointOutputDir);
    liblas::ReaderFactory rf;

    std::vector<std::string> files(config.pointFiles.begin(), config.pointFiles.end());
    FileSorter fs(10, -10);
    std::sort(files.begin(), files.end(), fs);
    std::list<Point_3> pts;
    std::list<liblas::Point> repeats;
    
    for (unsigned int i = 0; i < files.size(); ++i) {
        const std::string &pointFile = files[i];
        path oldPath(pointFile);
        path newPath(outDir / oldPath.filename());

        g_debug(" -- normalize (point) processing " << pointFile << " to " << newPath);
        if (exists(newPath)) {
            g_warn("The output file " << newPath << " already exists.");
            continue;
        }

        std::ifstream instr(pointFile.c_str(), std::ios::binary);
        liblas::Reader lasReader = rf.CreateWithStream(instr);
        liblas::Header lasHeader = lasReader.GetHeader();

        std::ofstream outstr(newPath.c_str(), std::ios::binary);
        liblas::Header outHeader(lasHeader);
        liblas::Writer lasWriter(outstr, outHeader);
        
        while (lasReader.ReadNextPoint()) {
            const liblas::Point &opt = lasReader.GetPoint();
            if(opt.GetClassification().GetClass() != 2)
                continue;
            Point_3 p(opt.GetX(), opt.GetY(), opt.GetZ());
            pts.push_back(std::move(p));
        }
        // TODO: This is not "correct".
        for(const liblas::Point &opt : repeats) {
            Point_3 p(opt.GetX(), opt.GetY(), opt.GetZ());
            pts.push_back(std::move(p));            
        }

        lasReader.Reset();
        
        Delaunay dt(pts.begin(), pts.end());
        pts.clear();
        repeats.clear();
        
        Face_handle hint;
        uint64_t ptCount = 0;
        uint64_t ptCountByReturn[255];
        for(int i = 0; i < 5; ++i)
            ptCountByReturn[i] = 0;
        
        while (lasReader.ReadNextPoint()) {
            const liblas::Point &opt = lasReader.GetPoint();
            if(opt.GetClassification().GetClass() == 2)
                continue;
            Point_3 p(opt.GetX(), opt.GetY(), opt.GetZ());
            hint = dt.locate(p, hint);
            if(dt.is_infinite(hint)) {
                repeats.push_back(opt);
                continue;
            }
            double area = 0.0;
            double total = 0.0;
            for(int i = 0; i < 3; ++i) {
                Point_3 p1 = hint->vertex(i)->point();
                Point_3 p2 = hint->vertex((i + 1) % 3)->point();
                Point_3 p3 = hint->vertex((i + 2) % 3)->point();
                double h = computeArea(p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), p.x(), p.y(), p.z());
                area += h;
                total += h * p3.z();
            }
            liblas::Point npt(opt);
            if(p.z() > 0.0 || !config.dropNegative)
                npt.SetZ(g_max(0.0, p.z() - total / area));
            lasWriter.WritePoint(npt);
            ptCountByReturn[npt.GetReturnNumber() - 1]++;
            ++ptCount;
        }

        g_debug(" -- setting original count " << outHeader.GetPointRecordsCount() << " to " << ptCount);
        outHeader.SetPointRecordsCount(ptCount);
        for(int i = 0; i < 5; ++i)
            outHeader.SetPointRecordsByReturnCount(i + 1, ptCountByReturn[i]);
        
        // Keep the bounds of the infinite faces to triangulate in the next round.
        for(All_faces_iterator it = dt.all_faces_begin(); it != dt.all_faces_end(); ++it) {
            if(dt.is_infinite(it)) {
                for(int i = 0; i < 3; ++i) {
                    if(it->vertex(i)->point().z() > 0)
                        pts.push_back(it->vertex(i)->point());
                }
            }
        }
        
        instr.close();
        outstr.close();
    }


}
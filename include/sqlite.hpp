#ifndef __SQLITE_HPP__
#define __SQLITE_HPP__

#include <sstream>
#include <iomanip>
#include <vector>

#include <sqlite3.h>
//#include <spatialite/gaiageo.h>
#include <spatialite.h>

#include "geotools.h"
#include "util.hpp"

namespace geotools {

    namespace db {

        class SQLite {
        private:
            int m_type;
            int m_srid;
            bool m_trans;
            std::string m_file;
            std::map<std::string, int> m_fields;

            sqlite3 *m_db;
            sqlite3_stmt *m_stmt;
            void *m_cache;

            bool exists(const std::string& name) {
                struct stat buffer;
                return (stat(name.c_str(), &buffer) == 0);
            }

        public:
            const static int INTEGER = 1;
            const static int DOUBLE = 2;
            const static int STRING = 3;
            const static int BLOB = 4;

            const static int POINT = 1;
            const static int LINESTRING = 2;
            const static int POLYGON = 3;

            SQLite(const std::string &file, int type, int srid,
                    const std::map<std::string, int> &fields) :
            m_type(type),
            m_srid(srid),
            m_trans(false),
            m_file(file),
            m_fields(fields),
            m_db(nullptr), m_stmt(nullptr), m_cache(nullptr) {
                init();
            }

            SQLite(const std::string &file, int type, int srid) :
            m_type(type),
            m_srid(srid),
            m_file(file),
            m_fields(std::map<std::string, int>()) {
                init();
            }

            SQLite(const std::string &file) :
            m_type(-1),
            m_srid(-1),
            m_file(file),
            m_fields(std::map<std::string, int>()) {
                init();
            }

            void clear();

            void addPoint(double x, double y, double z, const std::map<std::string, std::string> &fields) {

                std::stringstream ss;
                ss << std::setprecision(12) << "POINTZ(" << x << " " << y << " " << z << ")";
                std::string q = ss.str();

                sqlite3_reset(m_stmt);
                sqlite3_clear_bindings(m_stmt);
                sqlite3_bind_text(m_stmt, 1, q.c_str(), (int) q.size(), SQLITE_STATIC);
                int i = 1;
                for (auto it = fields.begin(); it != fields.end(); ++it) {
                    ++i;
                    switch (m_fields[it->first]) {
                        case SQLite::INTEGER:
                            sqlite3_bind_int(m_stmt, i, atoi(it->second.c_str()));
                            break;
                        case SQLite::DOUBLE:
                            sqlite3_bind_double(m_stmt, i, atof(it->second.c_str()));
                            break;
                        case SQLite::STRING:
                            sqlite3_bind_text(m_stmt, i, it->second.c_str(), (int) it->second.size(), SQLITE_STATIC);
                            break;
                    }
                }

                int ret = sqlite3_step(m_stmt);
                if (!(ret == SQLITE_DONE || ret == SQLITE_ROW))
                    handleError("Failed to insert row: ");
            }

            /**
             * Returns the number of points that can be added at one time by 
             * addPoints. This is computed from the parameter limit and the 
             * number of fields.
             */
            uint64_t maxAddPointCount() {
                uint64_t c = (uint64_t) (sqlite3_limit(m_db, SQLITE_LIMIT_VARIABLE_NUMBER, -1) / (m_fields.size() + 1));
                g_debug("max point count: " << c << ", " << sqlite3_limit(m_db, SQLITE_LIMIT_VARIABLE_NUMBER, -1));
                return c;
            }

            void makeFast() {
                g_debug("makeFast");
                char *err;
                if (SQLITE_OK != sqlite3_exec(m_db, "PRAGMA synchronous = OFF", NULL, NULL, &err))
                    handleError("Failed to set synchronous: ", err);
                if (SQLITE_OK != sqlite3_exec(m_db, "PRAGMA journal_mode = MEMORY", NULL, NULL, &err))
                    handleError("Failed to set journal mode: ", err);
                g_debug("done makeFast");
            }

            void makeSlow() {
                char *err;
                if (SQLITE_OK != sqlite3_exec(m_db, "PRAGMA synchronous = ON", NULL, NULL, &err))
                    handleError("Failed to set synchronous: ", err);
                if (SQLITE_OK != sqlite3_exec(m_db, "PRAGMA journal_mode = DELETE", NULL, NULL, &err))
                    handleError("Failed to set journal mode: ", err);
            }

            void setCacheSize(size_t size) {
                g_debug("setCacheSize " << size);
                char *err;
                std::stringstream ss;
                ss << "PRAGMA cache_size = " << size;
                if (SQLITE_OK != sqlite3_exec(m_db, ss.str().c_str(), NULL, NULL, &err))
                    handleError("Failed to set cache size: ", err);
                g_debug("done setCacheSize " << size);
            }

            void addPoints(std::vector<std::unique_ptr<geotools::util::Point> > &points) {
                for (const std::unique_ptr<geotools::util::Point> &pt : points)
                    addPoint(pt->x, pt->y, pt->z, pt->fields);
            }

            void dropGeomIndex() {
                char *err;
                if (SQLITE_OK != sqlite3_exec(m_db, "SELECT DisableSpatialIndex('data', 'geom')", NULL, NULL, &err))
                    handleError("Failed to drop geometry index: ", err);
            }

            void createGeomIndex() {
                char *err;
                if (SQLITE_OK != sqlite3_exec(m_db, "SELECT CreateSpatialIndex('data', 'geom'); VACUUM data", NULL, NULL, &err))
                    handleError("Failed to drop geometry index: ", err);
            }

            static int getPointsCallback(void *resultPtr, int cols, char **values, char **colnames) {
                using namespace geotools::util;
                std::vector<std::unique_ptr<Point> > *result = (std::vector<std::unique_ptr<Point> > *) resultPtr;
                Point *pt = new Point();
                for (int i = 0; i < cols; ++i) {
                    std::string colname(colnames[i]);
                    if (colname == "geomx") {
                        pt->x = atof(values[i]);
                    } else if (colname == "geomy") {
                        pt->y = atof(values[i]);
                    } else if (colname == "geomz") {
                        pt->z = atof(values[i]);
                    } else {
                        pt->fields[colname] = std::string(values[i]);
                    }
                }
                result->push_back(std::unique_ptr<Point>(pt));
                return 0;
            }

            void getPoints(std::vector<std::unique_ptr<geotools::util::Point> > &points,
                    const geotools::util::Bounds &bounds) {

                std::stringstream ss;
                ss << std::setprecision(12);
                ss << "SELECT X(geom) AS geomx, Y(geom) AS geomy, Z(geom) AS geomz";
                for (auto it = m_fields.begin(); it != m_fields.end(); ++it)
                    ss << ", " << it->first;
                ss << " FROM data WHERE Within(geom, GeomFromText('POLYGON((";
                ss << bounds.minx() << " " << bounds.miny() << ",";
                ss << bounds.maxx() << " " << bounds.miny() << ",";
                ss << bounds.maxx() << " " << bounds.maxy() << ",";
                ss << bounds.minx() << " " << bounds.maxy() << ",";
                ss << bounds.minx() << " " << bounds.miny();
                ss << "))', SRID(geom)))";

                std::string q = ss.str();
                g_debug("getPoints: " << q);
                char *err;
                begin();
                if (SQLITE_OK != sqlite3_exec(m_db, q.c_str(),
                        geotools::db::SQLite::getPointsCallback, &points, &err)) {
                    rollback();
                    handleError("Failed to execute query.", err);
                }
                rollback();
            }

            static int countCallback(void *resultPtr, int cols, char **values, char **colnames) {
                uint64_t *result = (uint64_t *) resultPtr;
                *result = atoi(values[0]);
                return 0;
            }

            void getGeomCount(uint64_t *count) {
                char *err;
                begin();
                if (SQLITE_OK != sqlite3_exec(m_db, "SELECT COUNT(*) FROM data",
                        SQLite::countCallback, count, &err)) {
                    rollback();
                    handleError("Failed to retrieve record count: ", err);
                }
                rollback();
            }

            void begin() {
                m_trans = true;
                char *err;
                if (SQLITE_OK != sqlite3_exec(m_db, "BEGIN TRANSACTION", NULL, NULL, &err))
                    handleError("Failed to start transaction: ", err);
            }

            void rollback() {
                if (m_trans) {
                    m_trans = false;
                    char *err;
                    if (SQLITE_OK != sqlite3_exec(m_db, "ROLLBACK TRANSACTION", NULL, NULL, &err))
                        handleError("Failed to start transaction: ", err);
                }
            }

            void commit() {
                if (m_trans) {
                    m_trans = false;
                    char *err;
                    if (SQLITE_OK != sqlite3_exec(m_db, "COMMIT TRANSACTION", NULL, NULL, &err))
                        handleError("Failed to commit transaction: ", err);
                }
            }

            ~SQLite() {
                commit();
                sqlite3_finalize(m_stmt);
                sqlite3_close(m_db);
            }

            void handleError(const std::string &msg, char *err = 0);
            void init();

        };

        void SQLite::handleError(const std::string &msg, char *err) {
            if (!err)
                err = (char *) sqlite3_errmsg(m_db);
            const std::string msg0(err);
            sqlite3_free(err);
            g_runerr(msg << msg0);
        }

        static int tableInfoCallback(void *resultPtr, int cols, char **values, char **colnames) {
            std::map<std::string, int> *fields = (std::map<std::string, int> *) resultPtr;
            std::pair<std::string, int> kv;
            for (int i = 0; i < cols; ++i) {
                std::string name(colnames[i]);
                if (name == "name") {
                    kv.first = std::string(values[i]);
                } else if (name == "type") {
                    const std::string type(values[i]);
                    if (type == "REAL") {
                        kv.second = SQLite::DOUBLE;
                    } else if (type == "TEXT") {
                        kv.second = SQLite::STRING;
                    } else if (type == "BLOB") {
                        kv.second = SQLite::BLOB;
                    } else if (type == "INTEGER") {
                        kv.second = SQLite::INTEGER;
                    } else {
                        kv.second = 0;
                    }
                }
            }
            if (kv.first != "geom" && kv.first != "gid") {
                g_debug("Field: " << kv.first << ", " << kv.second);
                fields->insert(kv);
            }
            return 0;
        }

        static int sridCallback(void *resultPtr, int cols, char **values, char **colnames) {
            int *srid = (int *) resultPtr;
            *srid = atoi(values[0]);
            return 0;
        }

        void SQLite::init() {
            bool dbExists = exists(m_file);
            g_debug("Initializing sqlite with exists: " << dbExists);

            char *err;

            g_debug("Opening...");
            if (SQLITE_OK != sqlite3_open_v2(m_file.c_str(), &m_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL))
                handleError("Failed to open DB: ");

            m_cache = spatialite_alloc_connection();
            spatialite_init_ex(m_db, m_cache, 0);

            if (dbExists) {
                g_debug("DB exists so collecting metadata.");
                if (SQLITE_OK != sqlite3_exec(m_db, "PRAGMA table_info('data')",
                        tableInfoCallback, &m_fields, &err))
                    handleError("Failed to read table info from database. Formatted incorrectly? ", err);
                g_debug(" -- get srid");
                if (SQLITE_OK != sqlite3_exec(m_db, "SELECT SRID(geom) AS geomsrid FROM data LIMIT 1",
                        sridCallback, &m_srid, &err))
                    handleError("Failed to read table info from database. Formatted incorrectly? ", err);
                g_debug("SRID: " << m_srid);
            } else {
                g_debug("Table is new, so initializing as spatial.");
                if (SQLITE_OK != sqlite3_exec(m_db, "SELECT InitSpatialMetadata(1)", NULL, NULL, &err))
                    handleError("Failed to initialize DB: ", err);
            }

            g_debug("Building DDL");

            std::stringstream ss;
            std::stringstream fn;
            std::stringstream fp;

            ss << "CREATE TABLE data(gid INTEGER PRIMARY KEY,";
            bool com = false;
            for (const auto &it : m_fields) {
                if (com) {
                    ss << ",";
                    fn << ",";
                    fp << ",";
                } else {
                    com = true;
                }

                fn << it.first;
                ss << it.first;
                fp << "?";

                switch (it.second) {
                    case SQLite::INTEGER:
                        ss << " integer";
                        break;
                    case SQLite::DOUBLE:
                        ss << " double";
                        break;
                    case SQLite::STRING:
                        ss << " text";
                        break;
                }
            }
            ss << ")";

            std::string q;

            if (!dbExists) {
                q.assign(ss.str());
                g_debug("Creating table: " << q);
                if (SQLITE_OK != sqlite3_exec(m_db, q.c_str(), NULL, NULL, &err))
                    handleError("Failed to create table: ", err);
                ss.str(std::string());
                ss.clear();
                ss << "SELECT AddGeometryColumn('data', 'geom', " << m_srid << ", '";
                switch (m_type) {
                    case POINT:
                        ss << "POINT";
                        break;
                    case LINESTRING:
                        ss << "LINESTRING";
                        break;
                    case POLYGON:
                        ss << "POLYGON";
                        break;
                }
                ss << "', 'XYZ')";
                q.assign(ss.str());

                g_debug("Creating geometry column: " << q);
                if (SQLITE_OK != sqlite3_exec(m_db, q.c_str(), NULL, NULL, &err))
                    handleError("Failed to create geometry: ", err);
            }

            ss.str(std::string());
            ss.clear();
            ss << "INSERT INTO data (geom, " << fn.str() << ") VALUES (GeomFromText(?, "
                    << m_srid << "), " << fp.str() << ")";

            q.assign(ss.str());
            g_debug("Preparing insert query: " << q);
            if (SQLITE_OK != sqlite3_prepare_v2(m_db, q.c_str(), (int) q.size(), &m_stmt, NULL))
                handleError("Failed to prepare insert statement: ");

        }

        void SQLite::clear() {
            g_debug("Deleting existing records.");
            char *err;
            begin();
            dropGeomIndex();
            if (SQLITE_OK != sqlite3_exec(m_db, "DELETE FROM data", NULL, NULL, &err)) {
                rollback();
                handleError("Failed to clear data table: ", err);
            }
            commit();
        }


    } // db

} // geotools

#endif


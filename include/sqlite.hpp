#ifndef __SQLITE_HPP__
#define __SQLITE_HPP__

#include <sstream>
#include <iomanip>
#include <vector>

#include <sqlite3.h>
#include <spatialite/gaiageo.h>
#include <spatialite.h>

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
				return (stat (name.c_str(), &buffer) == 0); 
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
				const std::map<std::string, int> &fields, bool clear = false) :
				m_file(file),
				m_type(type), 
				m_srid(srid),
				m_fields(fields) {
				init(clear);
			}

			SQLite(const std::string &file, int type, int srid, bool clear = false) :
				SQLite(file, type, srid, std::map<std::string, int>(), clear) {
			}
		
			SQLite(const std::string &file) :
				SQLite(file, -1, -1, std::map<std::string, int>(), false) {
			
			}
	
			void addPoint(double x, double y, double z, const std::map<std::string, std::string> &fields) {

				std::stringstream ss;
				ss << std::setprecision(12) << "POINTZ(" << x << " " << y << " " << z << ")";
				std::string q = ss.str();

				sqlite3_reset(m_stmt);
				sqlite3_clear_bindings(m_stmt);
				sqlite3_bind_text(m_stmt, 1, q.c_str(), q.size(), SQLITE_STATIC);
				int i = 1;
				for(auto it = fields.begin(); it != fields.end(); ++it) {
					++i;
					switch(m_fields[it->first]) {
					case SQLite::INTEGER:
						sqlite3_bind_int(m_stmt, i, atoi(it->second.c_str()));
						break;
					case SQLite::DOUBLE:
						sqlite3_bind_double(m_stmt, i, atof(it->second.c_str()));
						break;
					case SQLite::STRING:
						sqlite3_bind_text(m_stmt, i, it->second.c_str(), it->second.size(), SQLITE_STATIC);
						break;
					}
				}

				int ret = sqlite3_step(m_stmt);
				if(!(ret == SQLITE_DONE || ret == SQLITE_ROW))
					handleError("Failed to insert row: ");
			}

			static int getPointsCallback(void *resultPtr, int cols, char **values, char **colnames) {
				using namespace geotools::util;
				std::vector<std::unique_ptr<Point> > *result = (std::vector<std::unique_ptr<Point> > *) resultPtr;
				Point *pt = new Point();
				for(int i = 0; i < cols; ++i) {
					std::string colname(colnames[i]);
					if(colname == "geomx") {
						pt->x = atof(values[i]);
					} else if(colname == "geomy") {
						pt->y = atof(values[i]);
					} else if(colname == "geomz") {
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
				
				begin();
				std::stringstream ss;
				ss << std::setprecision(12);
				ss << "SELECT X(geom) AS geomx, Y(geom) AS geomy, Z(geom) AS geomz";
				for(auto it = m_fields.begin(); it != m_fields.end(); ++it)
					ss << ", " << it->first;
				ss << " FROM data WHERE Within(geom, GeomFromText('POLYGON((";
				ss << bounds.minx() << " " << bounds.miny() << ",";
				ss << bounds.maxx() << " " << bounds.miny() << ",";
				ss << bounds.maxx() << " " << bounds.maxy() << ",";
				ss << bounds.minx() << " " << bounds.maxy() << ",";
				ss << bounds.minx() << " " << bounds.miny();
				ss << "))', SRID(geom)))";

				std::string q = ss.str();
				g_trace("getPoints: " << q);
				char *err;	
				if(SQLITE_OK != sqlite3_exec(m_db, q.c_str(),  
					geotools::db::SQLite::getPointsCallback, &points, &err))
					handleError("Failed to execute query.", err);
				rollback();
			}

			void begin() {
				m_trans = true;
				char *err;
				if(SQLITE_OK != sqlite3_exec(m_db, "BEGIN TRANSACTION", NULL, NULL, &err))
					handleError("Failed to start transaction: ", err);
			}

			void rollback() {
				if(m_trans) {
					m_trans = false;
					char *err;
					if(SQLITE_OK != sqlite3_exec(m_db, "ROLLBACK TRANSACTION", NULL, NULL, &err))
						handleError("Failed to start transaction: ", err);
				}
			}

			void commit() {
				if(m_trans) {
					m_trans = false;
					char *err;
					if(SQLITE_OK != sqlite3_exec(m_db, "COMMIT TRANSACTION", NULL, NULL, &err))
						handleError("Failed to commit transaction: ", err);
				}
			}

			~SQLite() {
				commit();
				sqlite3_finalize(m_stmt);
				sqlite3_close(m_db);
			}

			void handleError(const std::string &msg, char *err = 0);
			void init(bool clear = false);

		};

		void SQLite::handleError(const std::string &msg, char *err) {
			if(!err)
				err = (char *) sqlite3_errmsg(m_db);
			const std::string msg0(err);
			sqlite3_free(err);
			g_runerr(msg << msg0);
		}

		static int strToType(const std::string &type) {

		}
			
		static int tableInfoCallback(void *resultPtr, int cols, char **values, char **colnames) {
			std::map<std::string, int> *fields = (std::map<std::string, int> *) resultPtr;
			std::pair<std::string, int> kv;
			for(int i = 0; i < cols; ++i) {
				std::string name(colnames[i]);
				if(name == "name") {
					kv.first = std::string(values[i]);
				} else if(name == "type") {
					const std::string type(values[i]);
					if(type == "REAL") {
						kv.second = SQLite::DOUBLE;
					} else if(type == "TEXT") {
						kv.second = SQLite::STRING;
					} else if(type == "BLOB") {
						kv.second = SQLite::BLOB;
					} else if(type == "INTEGER") {
						kv.second = SQLite::INTEGER;
					} else {
						kv.second = 0;
					}
				}
			}
			if(kv.first != "geom") {
				g_trace("Field: " << kv.first << ", " << kv.second);
				fields->insert(kv);
			}
			return 0;
		}

		static int sridCallback(void *resultPtr, int cols, char **values, char **colnames) {
			int *srid = (int *) resultPtr;
			*srid = atoi(values[0]);
			return 0;
		}

		void SQLite::init(bool clear) {
			bool dbExists = exists(m_file);
			g_trace("Initializing sqlite with clear: " << clear << ", and exists: " << dbExists);
			
			char *err;

			g_trace("Opening...");
			if(SQLITE_OK != sqlite3_open_v2(m_file.c_str(), &m_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL))
				handleError("Failed to open DB: ");
			
			m_cache = spatialite_alloc_connection();
			spatialite_init_ex(m_db, m_cache, 0);

			if(dbExists) {
				g_trace("DB exists so collecting metadata.");
				if(SQLITE_OK != sqlite3_exec(m_db, "PRAGMA table_info('data')", 
					tableInfoCallback, &m_fields, &err))
					handleError("Failed to read table info from database. Formatted incorrectly? ", err);
				if(SQLITE_OK != sqlite3_exec(m_db, "SELECT SRID(geom) AS geomsrid FROM data", 
					sridCallback, &m_srid, &err))
					handleError("Failed to read table info from database. Formatted incorrectly? ", err);
				g_trace("SRID: " << m_srid);
			} else {
				g_trace("Table is new, so initializing as spatial.");
				if(SQLITE_OK != sqlite3_exec(m_db, "SELECT InitSpatialMetadata(1)", NULL, NULL, &err))
					handleError("Failed to initialize DB: ", err);
			}


			std::stringstream ss;
			std::stringstream fn;
			std::stringstream fp;

			ss << "CREATE TABLE data(";
			bool com = false;
			for(auto it = m_fields.begin(); it != m_fields.end(); ++it) {
				if(com) {
					ss << ", ";
					fn << ", ";
					fp << ", ";
				} else {
					com = true;
				}

				fn << it->first;
				ss << it->first;
				fp << "?";

				switch(it->second) {
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

			if(!dbExists) {
				q.assign(ss.str());
				g_trace("Creating table: " << q);
				if(SQLITE_OK != sqlite3_exec(m_db, q.c_str(), NULL, NULL, &err)) 
					handleError("Failed to create table: ", err);
				ss.str(std::string());
				ss.clear();
				ss << "SELECT AddGeometryColumn('data', 'geom', " << m_srid << ", '";
				switch(m_type) {
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

				g_trace("Creating geometry column: " << q);
				if(SQLITE_OK != sqlite3_exec(m_db, q.c_str(), NULL, NULL, &err))
					handleError("Failed to create geometry: ", err);
			}

			if(clear) {
				g_trace("Deleting existing records.");
				if(SQLITE_OK != sqlite3_exec(m_db, "DELETE FROM data", NULL, NULL, &err))
					handleError("Failed to clear data table: ", err);
			}

			ss.str(std::string());
			ss.clear();
			ss << "INSERT INTO data (geom, " << fn.str() << ") VALUES (GeomFromText(?, " 
				<< m_srid << "), " << fp.str() << ")";

			q.assign(ss.str());
			g_trace("Preparing insert query: " << q);
			if(SQLITE_OK != sqlite3_prepare_v2(m_db, q.c_str(), q.size(), &m_stmt, NULL))
				handleError("Failed to prepare insert statement: ");

		}

	} // db

} // geotools

#endif


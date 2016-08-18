#ifndef __SQLITE_HPP__
#define __SQLITE_HPP__

#include <sqlite3.h>
#include <spatialite/gaiageo.h>
#include <spatialite.h>

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
			
			const static int POINT = 1;
			const static int LINESTRING = 2;
			const static int POLYGON = 3;

			SQLite(const std::string &file, int type, int srid, const std::map<std::string, int> &fields, bool clear = false) :
				m_file(file),
				m_type(type), 
				m_srid(srid),
				m_fields(fields) {
				init(clear);
			}

			SQLite(const std::string &file, int type, int srid, bool clear = false) :
				SQLite(file, type, srid, std::map<std::string, int>(), clear) {
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

			void begin() {
				m_trans = true;
				char *err;
				if(SQLITE_OK != sqlite3_exec(m_db, "BEGIN TRANSACTION", NULL, NULL, &err))
					handleError("Failed to start transaction: ", err);
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

		void SQLite::init(bool clear) {

			bool doInit = !exists(m_file);

			char *err;

			if(SQLITE_OK != sqlite3_open_v2(m_file.c_str(), &m_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL))
				handleError("Failed to open DB: ");
			
			m_cache = spatialite_alloc_connection();
			spatialite_init_ex(m_db, m_cache, 0);

			if(doInit) {
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
					ss << ",";
					fn << ",";
					fp << ",";
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

			if(doInit) {
				q.assign(ss.str());
				if(SQLITE_OK != sqlite3_exec(m_db, q.c_str(), NULL, NULL, &err)) 
					handleError("Failed to create table: ", err);

				ss.str(std::string());
				ss.clear();
				ss << "SELECT AddGeometryColumn('data', 'Geometry', " << m_srid << ", '";
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
				if(SQLITE_OK != sqlite3_exec(m_db, q.c_str(), NULL, NULL, &err))
					handleError("Failed to create geometry: ", err);

			}

			if(clear) {
				if(SQLITE_OK != sqlite3_exec(m_db, "DELETE FROM data", NULL, NULL, &err))
					handleError("Failed to clear data table: ", err);
			}

			ss.str(std::string());
			ss.clear();
			ss << "INSERT INTO data (Geometry, " << fn.str() << ") VALUES (GeomFromText(?, " << m_srid << "), " << fp.str() << ")";

			q.assign(ss.str());
			if(SQLITE_OK != sqlite3_prepare_v2(m_db, q.c_str(), q.size(), &m_stmt, NULL))
				handleError("Failed to prepare insert statement: ");

		}

	} // db

} // geotools

#endif


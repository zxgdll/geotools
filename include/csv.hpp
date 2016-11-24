/* 
 * File:   csv.hpp
 * Author: robskelly
 *
 * Created on November 24, 2016, 3:17 PM
 */

#ifndef CSV_HPP
#define	CSV_HPP

#include <string>
#include <memory>
#include <vector>
#include <fstream>
#include <sstream>

#ifdef	__cplusplus
extern "C" {
#endif

namespace geotools {
    namespace csv {

        class CSVReader {
        private:
            std::string m_filename;
            std::unique_ptr<std::ifstream> m_str;
            std::vector<std::string> m_header;
            char *m_buf;
            uint32_t m_buflen;
            uint8_t m_type;

            void parseBuf(std::vector<std::string> &out) {
                out.clear();
                if(m_buf[0] == '#')
                    return;
                bool esc = false;
                bool quote = false;
                int i = 0;
                char c;
                std::stringstream ss;
                while((c = m_buf[i++]) != '\0') {
                    switch(c) {
                        case '\\':
                            if(!esc) {
                                esc = true;
                            } else {
                                ss << c;
                                esc = false;
                            }
                            break;
                        case '"':
                            if(!esc) {
                                quote = !quote;
                            } else {
                                ss << c;
                                esc = false;
                            }
                            break;
                        case ',':
                            out.push_back(ss.str());
                            ss.clear();
                            ss.str(std::string());
                            break;
                        default:
                            ss << c;
                            break;
                    }
                }
                out.push_back(ss.str());
            }
        public:
            CSVReader(const std::string &filename, uint32_t buflen = 2048) : 
                m_filename(filename),
                m_buf(nullptr),
                m_buflen(buflen) {
            }
            bool next(std::unordered_map<std::string, std::string> &row) {
                if(!m_str.get()) {
                    m_buf = (char *) std::malloc(m_buflen + 1);
                    g_debug(m_filename);
                    m_str.reset(new std::ifstream(m_filename));
                    m_str->getline(m_buf, m_buflen);
                    g_debug(" -- fail: " << m_str->fail() << "; " << m_str->eof());
                    if(m_str->fail() || m_str->eof())
                        g_runerr("Failed to parse file: " << m_filename);
                    parseBuf(m_header);
                }
                m_str->getline(m_buf, m_buflen);
                if(m_str->fail() || m_str->eof())
                        return false;
                std::vector<std::string> values;
                parseBuf(values);
                if(values.size() < m_header.size())
                    return false;
                for(uint32_t i = 0; i < m_header.size(); ++i)
                    row[m_header[i]] = values[i];
                return true;
            }
            ~CSVReader() {
                if(m_buf)
                    delete m_buf;
            }
        };

    }
}

#ifdef	__cplusplus
}
#endif

#endif	/* CSV_HPP */


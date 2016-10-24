#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <memory>
#include <cstdio>
#include <cerrno>

#include "omp.h"

#include "liblas/liblas.hpp"

#include "geotools.h"
#include "sortedpointstream.hpp"
#include "lasutil.hpp"

using namespace geotools::las;
using namespace geotools::util;

double LASPoint::scaleX = 0;
double LASPoint::scaleZ = 0;
double LASPoint::scaleY = 0;

std::string _rowFile(unsigned int row) {
	std::stringstream fs;
	fs << "/tmp/sps_" << row << ".tmp";
	return fs.str();
}

unsigned int _getRow(double y, const Bounds &bounds, double blockSize) {
	if(blockSize < 0) {
		return (unsigned int) ((y - bounds.maxy()) / blockSize);
	} else {
		return (unsigned int) ((bounds.miny() - y) / blockSize);
	}
}

SortedPointStream::SortedPointStream(const std::list<std::string> &files, double blockSize, bool rebuild) {
	init(files, blockSize, rebuild);
}

SortedPointStream::~SortedPointStream() {
	delete m_pf;
	/*
	for(const unsigned int &row : m_rows) {
		try {
			std::remove(_rowFile(row).c_str());
		} catch(...) {
			g_warn("Failed to delete " << _rowFile(row));
		}
	}
	*/
}

unsigned int SortedPointStream::pointCount() const {
	return m_pointCount;
}

unsigned int SortedPointStream::rowCount() const {
	return m_rows.size();
}

void SortedPointStream::init(const std::list<std::string> &files, double blockSize, bool rebuild) {

	m_bounds.collapse();
	m_pointCount = 0;
	m_row = 0;

	liblas::ReaderFactory rf;
	std::unordered_map<unsigned int, std::string> rowFileNames;

	for(const std::string &file : files) {
		g_debug(" -- init: opening file: " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		LASPoint::setScale(lasHeader.GetScaleX(), lasHeader.GetScaleY(), lasHeader.GetScaleZ());

		g_debug(" -- init computing bounds");
		Bounds fileBounds;
		if(!LasUtil::computeLasBounds(lasHeader, fileBounds, 2))
			LasUtil::computeLasBounds(lasReader, fileBounds, 2); // If the header bounds are bogus.
		m_bounds.extend(fileBounds);
	}
	m_bounds.snap(g_abs(blockSize));

	m_pf = new PointFile("/tmp/pf.tmp", m_bounds, blockSize);
	m_pf->openWrite();

	for(const std::string &file : files) {
		g_debug(" -- init: opening file: " << file);
		std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
		liblas::Reader lasReader = rf.CreateWithStream(instr);
		liblas::Header lasHeader = lasReader.GetHeader();

		LASPoint pt;
		while(lasReader.ReadNextPoint()) {
			pt.update(lasReader.GetPoint());
			m_pf->addPoint(pt);
		}
	}

	m_pf->close();
	m_pf->openRead();

		/*
	if(rebuild) {
	
		std::vector<std::string> files0(files.begin(), files.end());
		std::unordered_map<unsigned int, unsigned int> rowCounts;
		std::unordered_map<unsigned int, std::FILE*> rowFiles;
		std::unordered_map<unsigned int, unsigned long> rowTimes; // row --> time
		std::map<unsigned long, unsigned int> timeRows; // time --> row

		unsigned int pointCount = 0;
		unsigned long openFileId = 0;

		#pragma omp parallel reduction(+:pointCount)
		{

			#pragma omp for
			for(unsigned int i = 0; i < files0.size(); ++i) {
				
				const std::string &file = files0[i];

				g_debug(" -- init: opening file: " << file << "; " << omp_get_thread_num());
				std::ifstream instr(file.c_str(), std::ios::in|std::ios::binary);
				liblas::Reader lasReader = rf.CreateWithStream(instr);
				liblas::Header lasHeader = lasReader.GetHeader();

				g_debug(" -- init writing points" << "; " << omp_get_thread_num());
				LASPoint pt;
				while(lasReader.ReadNextPoint()) {
					pt.update(lasReader.GetPoint());
					unsigned int row = _getRow(pt.y, m_bounds, blockSize);

					if(rowFiles.find(row) == rowFiles.end()) {
						#pragma omp critical(__sps_update)
						{
							if(rowFiles.find(row) == rowFiles.end()) {
								++openFileId;
								std::string rowFile = _rowFile(row);
								FILE *os;
								while(!(os = std::fopen(rowFile.c_str(), "wb"))) {
									g_debug(" -- closing file -- limit reached");
									if(!rowFiles.size())
										g_runerr("Closed all files but couldn't open a new one... something's wrong.");
									const auto &rt = timeRows.begin(); // Get the oldest item.
									if(std::fclose(rowFiles[rt->second]))
										g_warn("Failed to close file.");
									rowFiles.erase(rt->second);
									timeRows.erase(rt->first);
									sleep(1);
								}
								std::fseek(os, sizeof(unsigned int), SEEK_SET);
								rowFiles[row] = os;
								timeRows[openFileId] = row;
								rowTimes[row] = openFileId;
								rowCounts[row] = 0;
								m_rows.push_back(row);
							}
						}

						pt.write(rowFiles[row]);
					}


					#pragma omp critical(__sps_update)
					{
						openFileId++;
						timeRows.erase(rowTimes[row]);
						rowTimes[row] = openFileId;
						timeRows[openFileId] = row;
						rowCounts[row]++;
					}

					++pointCount;
				}
			}

			m_pointCount = pointCount;
		}

		for(const auto &it : rowFiles) {
			unsigned int c = rowCounts[it.first];
			// it.second->seekp(0, std::ios::beg);
			std::fseek(it.second, 0, SEEK_SET);
			//it.second->write((char *) &c, sizeof(unsigned int));
			std::fwrite((const void *) &c, sizeof(unsigned int), 1, it.second);
			std::fclose(it.second);
		}

	} // rebuild
	*/
}

const Bounds& SortedPointStream::bounds() const {
	return m_bounds;
}

bool SortedPointStream::next(std::list<std::shared_ptr<LASPoint> > &pts) {
	return m_pf->nextRow(pts);

	/*
	#pragma omp critical(__sps_init)
	{
		if(m_row < m_rows.size()) {

			unsigned int row = m_rows[m_row];

			std::ifstream instr(_rowFile(row).c_str(), std::ios::binary);
			unsigned int count;
			instr.read((char *) &count, sizeof(unsigned int));
			if(count) {
				for(unsigned int i = 0; i < count; ++i) {
					std::shared_ptr<LASPoint> pt(new LASPoint());
					pt->read(instr);
					pts.push_back(std::move(pt));
				}
			}

			++m_row;
		}
	}
	*/
}

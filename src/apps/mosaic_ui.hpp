#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>
#include <vector>
#include <string>
#include <cmath>

#include <QtWidgets/QWidget>
#include <QDir>
#include <QTCore>
#include <QMessageBox>

#include "mosaic.hpp"
#include "util.hpp"
#include "ui_mosaic.h"

namespace geotools {

	namespace raster {

		class MosaicCallbacks : public QObject, public geotools::util::Callbacks {
			Q_OBJECT
		public:
			void fileCallback(float status) {
				emit fileProgress((int) std::round(status * 100));
			}

			void overallCallback(float status) {
				emit overallProgress((int) std::round(status * 100));
			}

		signals:
			void fileProgress(int);
			void overallProgress(int);
		};

	}

	namespace ui {

		class WorkerThread : public QThread {
		public:
			void init(QWidget *parent, geotools::raster::MosaicCallbacks *callbacks, const std::vector<std::string> &tifFiles, 
				const std::string &destFile, double distance, int tileSize, int threads) {
				m_callbacks = callbacks;
				m_tifFiles = tifFiles;
				m_destFile = destFile;
				m_distance = distance;
				m_tileSize = tileSize;
				m_threads = threads;
				m_parent = parent;
			}

		private:
			std::vector<std::string> m_tifFiles;
			std::string m_destFile;
			double m_distance;
			int m_tileSize;
			int m_threads;
			geotools::raster::MosaicCallbacks *m_callbacks;
			QWidget *m_parent;

			void run() {
				geotools::raster::Mosaic m;
				m.setCallbacks(m_callbacks);
				try {
					m.mosaic(m_tifFiles, m_destFile, m_distance, m_tileSize, m_threads);
				} catch(const std::exception &e) {
					QMessageBox err(m_parent);
					err.setText("Error");
					err.setInformativeText(QString(e.what()));
					err.exec();
				}
			}
		};

		class MosaicForm : public QWidget, public Ui::MosaicForm {
			Q_OBJECT
		private:
			QWidget *m_form;
			std::string m_destFile;
			std::vector<std::string> m_tifFiles;
			QDir *m_last;
			int m_tileSize;
			double m_distance;
			int m_threads;
			void updateFileList();
			void updateFileButtons();
			void checkRun();
			geotools::raster::MosaicCallbacks *m_callbacks;
			WorkerThread m_workerThread;

		public:
			MosaicForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *Form);
			~MosaicForm();

		public slots:
			void fileListSelectionChanged();
			void selectFilesClicked();
			void removeFilesClicked();
			void clearFilesClicked();
			void cancelClicked();
			void runClicked();
			void destFileClicked();
			void upClicked();
			void downClicked();
			void distanceChanged(QString);
			void threadsChanged(int);
			void tileSizeChanged(QString);
			void done();
		};

	}

}

#endif


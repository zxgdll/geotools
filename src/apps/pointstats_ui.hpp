#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>
#include <cmath>

#include <QtWidgets/QWidget>
#include <QDir>
#include <QMessageBox>
#include <QtCore>

#include "util.hpp"
#include "pointstats.hpp"
#include "ui_pointstats.h"

namespace geotools {

	namespace point {


		class PointStatsCallbacks : public QObject, public geotools::util::Callbacks {
			Q_OBJECT
		public:
			void fileCallback(float status) const {
				emit fileProgress((int) std::round(status * 100));
			}

			void overallCallback(float status) const {
				emit overallProgress((int) std::round(status * 100));
			}

		signals:
			void fileProgress(int) const;
			void overallProgress(int) const;
		};

	}

	namespace ui {


		class WorkerThread;

		class PointStatsForm : public QWidget, public Ui::PointStatsForm {
		friend class WorkerThread;
			Q_OBJECT
		private:
			QWidget *m_form;
			std::string m_destFile;
			std::list<std::string> m_sourceFiles;
			std::set<unsigned char> m_classes;
			bool m_fill;
			bool m_snap;
			int m_vsrid;
			int m_hsrid;
			int m_attribute;
			int m_type;
			int m_quantile;
			int m_quantiles;
			double m_resolution;
			unsigned int m_threads;
			unsigned char m_angleLimit;
			unsigned char m_gapFunction;
			QDir m_last;
			geotools::util::Callbacks *m_callbacks;
			WorkerThread *m_workerThread;

			void updateFileList();
			void updateFileButtons();
			void updateTypeUi();
			void checkRun();
			
		public:
			PointStatsForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *Form);
			~PointStatsForm();

		public slots:
			void fileListSelectionChanged();
			void selectFilesClicked();
			void removeFilesClicked();
			void clearFilesClicked();
			void cancelClicked();
			void runClicked();
			void destFileClicked();
			void crsConfigClicked();
			void typeSelected(int);
			void threadsChanged(int);
			void quantileChanged(int);
			void quantilesChanged(int);
			void attributeSelected(int);
			void snapToGridChanged(bool);
			void resolutionChanged(double);
			void maxAngleChanged(int);
			void gapFunctionSelected(int);
			void done();

		};


		class WorkerThread : public QThread {
		private:
			geotools::util::Bounds m_bounds;
			PointStatsForm *m_parent;

			void run() {
				try {
					geotools::point::PointStats l;
					
					geotools::point::PointStatsConfig config;
					config.dstFile = m_parent->m_destFile;
					config.sourceFiles = m_parent->m_sourceFiles;
					config.classes = m_parent->m_classes;
					config.hsrid = m_parent->m_hsrid;
					config.attribute = m_parent->m_attribute;
					config.type = m_parent->m_type;
					config.resolution = m_parent->m_resolution;
					config.bounds = m_bounds;
					config.angleLimit = m_parent->m_angleLimit;
					config.fill = m_parent->m_fill;
					config.snap = m_parent->m_snap;
					config.threads = m_parent->m_threads;
					config.gapFractionType = m_parent->m_gapFunction;

					l.pointstats(config, m_parent->m_callbacks);
				} catch(const std::exception &e) {
					QMessageBox err((QWidget *) m_parent);
					err.setText("Error");
					err.setInformativeText(QString(e.what()));
					err.exec();
				}
			}
		public:
			void init(PointStatsForm *parent, const geotools::util::Bounds &bounds) {
				m_parent = parent;
				m_bounds = bounds;
			}
		};

	}

}

#endif


#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>
#include <cmath>

#include <QtWidgets/QWidget>
#include <QDir>
#include <QMessageBox>
#include <QTCore>

#include "util.hpp"
#include "lasgrid.hpp"
#include "ui_lasgrid.h"

namespace geotools {

	namespace las {


		class LasCallbacks : public QObject, public geotools::util::Callbacks {
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
			void init(QWidget *parent, geotools::las::LasCallbacks *callbacks, const std::vector<std::string> &lasFiles, 
				const std::string &destFile, const std::set<int> &classes, int hsrid, int attribute, int type, 
				double radius, double resolution, geotools::util::Bounds &bounds, unsigned char angleLimit, bool fill, bool snap) {
				m_callbacks = callbacks;
				m_lasFiles = lasFiles;
				m_destFile = destFile;
				m_classes = classes;
				m_hsrid = hsrid;
				m_attribute = attribute;
				m_type = type;
				m_radius = radius;
				m_resolution = resolution;
				m_bounds = bounds;
				m_angleLimit = angleLimit;
				m_fill = fill;
				m_snap = snap;
			}

		private:
			std::vector<std::string> m_lasFiles;
			std::string m_destFile;
			std::set<int> m_classes;
			int m_hsrid;
			int m_attribute;
			int m_type;
			double m_radius;
			double m_resolution;
			geotools::util::Bounds m_bounds;
			unsigned char m_angleLimit;
			bool m_fill;
			bool m_snap;
			geotools::las::LasCallbacks *m_callbacks;
			QWidget *m_parent;

			void run() {
				geotools::las::LasGrid l;
				l.setCallbacks(m_callbacks);
				try {
					l.lasgrid(m_destFile, m_lasFiles, m_classes, m_hsrid, m_attribute, m_type, m_radius, m_resolution, 
						m_bounds, m_angleLimit, m_fill, m_snap);
				} catch(const std::exception &e) {
					QMessageBox err(m_parent);
					err.setText("Error");
					err.setInformativeText(QString(e.what()));
					err.exec();
				}
			}
		};


		class LasgridForm : public QWidget, public Ui::LasgridForm {
			Q_OBJECT
		private:
			QWidget *m_form;
			std::string m_destFile;
			std::vector<std::string> m_lasFiles;
			std::set<int> m_classes;
			int m_vsrid;
			int m_hsrid;
			int m_attribute;
			int m_type;
			unsigned char m_angleLimit;
			bool m_fill;
			bool m_snap;
			double m_resolution;
			double m_radius;
			int m_quantile;
			int m_quantiles;
			QDir *m_last;
			geotools::las::LasCallbacks *m_callbacks;
			WorkerThread m_workerThread;

			void updateFileList();
			void updateFileButtons();
			void checkRun();

		public:
			LasgridForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *Form);
			~LasgridForm();

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
			void quantileChanged(int);
			void quantilesChanged(int);
			void attributeSelected(int);
			void radiusChanged(double);
			void snapToGridChanged(bool);
			void resolutionChanged(double);
			void maxAngleChanged(int);

		};

	}

}

#endif


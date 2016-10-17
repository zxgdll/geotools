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
			unsigned int m_quantileFilter;
			unsigned int m_quantileFilterFrom;
			unsigned int m_quantileFilterTo;
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

			void destFileClicked();
			void snapToGridChanged(bool);
			
			void cancelClicked();
			void runClicked();

			void crsConfigClicked();

			void typeSelected(int);

			void threadsChanged(int);

			void quantileChanged(int);
			void quantilesChanged(int);

			void attributeSelected(int);
			void resolutionChanged(double);
			void gapFunctionSelected(int);

			void quantileFilterFromChanged(int);
			void quantileFilterToChanged(int);
			void quantileFilterChanged(int);
			void maxAngleChanged(int);
			void classItemClicked(QListWidgetItem*);

			void done();

		};


		class WorkerThread : public QThread {
		private:
			geotools::util::Bounds m_bounds;
			PointStatsForm *m_parent;
			std::string m_error;

			void run();
		public:
			void init(PointStatsForm *parent, const geotools::util::Bounds &bounds) {
				m_parent = parent;
				m_bounds = bounds;
			}

			bool hasError();

			std::string getError();

		};

	}

}

#endif


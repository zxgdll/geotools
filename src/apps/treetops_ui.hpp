#ifndef __TREETOPS_UI_HPP__
#define __TREETOPS_UI_HPP__

#include <QtWidgets/QWidget>
#include <QMessageBox>
#include <QtCore>
#include <QDir>

#include "util.hpp"
#include "treetops.hpp"
#include "ui_treetops.h"

using namespace geotools::treetops::config;

namespace geotools {

	namespace trees {


		class TreetopsCallbacks : public QObject, public geotools::util::Callbacks {
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

		class WorkerThread;

		class TreetopsForm : public QWidget, public Ui::TreetopsForm {
		friend class WorkerThread;
			Q_OBJECT
		private:
			QWidget *m_form;
			QDir m_last;
			TreetopsConfig m_config;
			geotools::util::Callbacks *m_callbacks;
			WorkerThread *m_workerThread;

			void checkRun();

		public:
			TreetopsForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *Form);
			~TreetopsForm();

		public slots:
			void doSmoothChanged(bool);
			void doTopsChanged(bool);
			void doCrownsChanged(bool);

			void smoothWindowSizeChanged(int);
			void smoothStdDevChanged(double);
			void smoothInputFileChanged(QString);
			void smoothOutputFileChanged(QString);

			void topsMinHeightChanged(double);
			void topsWindowSizeChanged(int);
			void topsInputFileChanged(QString);
			void topsOutputFileChanged(QString);

			void crownsRadiusChanged(double);
			void crownsHeightFractionChanged(double);
			void crownsMinHeightChanged(double);
			void crownsTreetopsFileChanged(QString);
			void crownsInputFileChanged(QString);
			void crownsOutputRasterChanged(QString);
			void crownsOutputVectorChanged(QString);

			void exitClicked();
			void runClicked();
			void cancelClicked();

			void smoothInputFileClicked();
			void smoothOutputFileClicked();

			void topsInputFileClicked();
			void topsOutputFileClicked();
			
			void crownsOutputVectorClicked();
			void crownsOutputRasterClicked();
			void crownsTreetopsFileClicked();
			void crownsInputFileClicked();

			void done();
		};


		class WorkerThread : public QThread {
		private:
			TreetopsForm *m_parent;
			void run() {
				geotools::treetops::Treetops t;
				t.setCallbacks(m_parent->m_callbacks);
				try {

				} catch(const std::exception &e) {
					QMessageBox err((QWidget *) m_parent);
					err.setText("Error");
					err.setInformativeText(QString(e.what()));
					err.exec();
				}
			}
		public:
			void init(TreetopsForm *parent) {
				m_parent = parent;
			}
		};

	}

}

#endif


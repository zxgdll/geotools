#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <set>
#include <vector>
#include <string>
#include <math.h>

#include <QtWidgets/QWidget>
#include <QDir>
#include <QtCore>
#include <QMessageBox>

#include "mosaic.hpp"
#include "util.hpp"
#include "ui_mosaic.h"

#ifdef _MSC_VER
namespace std {

    inline double round(double value) {
        return value < 0 ? -std::floor(0.5 - value) : std::floor(0.5 + value);
    }
}
#endif

namespace geotools {

    namespace raster {

        class MosaicCallbacks : public QObject, public geotools::util::Callbacks {

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

        class MosaicForm : public QWidget, public Ui::MosaicForm {
            Q_OBJECT
            friend class WorkerThread;
        private:
            QWidget *m_form;
            std::string m_destFile;
            std::vector<std::string> m_tifFiles;
            QDir m_last;
            int m_tileSize;
            double m_distance;
            int m_threads;
            void updateFileList();
            void updateFileButtons();
            void checkRun();
            geotools::raster::MosaicCallbacks *m_callbacks;
            WorkerThread *m_workerThread;

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

        class WorkerThread : public QThread {
        private:
            MosaicForm *m_parent;

            void run() {
                geotools::raster::Mosaic m;
                m.setCallbacks(m_parent->m_callbacks);
                try {
                    m.mosaic(m_parent->m_tifFiles, m_parent->m_destFile, m_parent->m_distance,
                            m_parent->m_tileSize, m_parent->m_threads);
                } catch (const std::exception &e) {
                    QMessageBox err((QWidget *) m_parent);
                    err.setText("Error");
                    err.setInformativeText(QString(e.what()));
                    err.exec();
                }
            }

        public:

            void init(MosaicForm *parent) {
                m_parent = parent;
            }
        };


    }

}

#endif


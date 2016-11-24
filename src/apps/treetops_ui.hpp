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

    namespace treetops {

        class TreetopsCallbacks : public QObject, public geotools::util::Callbacks {
            Q_OBJECT
        public:
            TreetopsCallbacks() : geotools::util::Callbacks() {
            }
            void stepCallback(float status) const {
                emit stepProgress((int) std::round(status * 100));
            }

            void overallCallback(float status) const {
                emit overallProgress((int) std::round(status * 100));
            }

        signals:
            void stepProgress(int) const;
            void overallProgress(int) const;
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
            void updateView();

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
            void smoothOriginalCHMChanged(QString);
            void smoothSmoothedCHMChanged(QString);

            void topsMinHeightChanged(double);
            void topsWindowSizeChanged(int);
            void topsOriginalCHMChanged(QString);
            void topsSmoothedCHMChanged(QString);
            void topsTreetopsDatabaseChanged(QString);
            void topsTreetopsSRIDClicked();
            void topsTreetopsSRIDChanged(int);

            void crownsRadiusChanged(double);
            void crownsHeightFractionChanged(double);
            void crownsMinHeightChanged(double);
            void crownsTreetopsDatabaseChanged(QString);
            void crownsSmoothedCHMChanged(QString);
            void crownsCrownsRasterChanged(QString);
            void crownsCrownsDatabaseChanged(QString);

            void exitClicked();
            void runClicked();
            void cancelClicked();
            void helpClicked();
            
            void smoothOriginalCHMClicked();
            void smoothSmoothedCHMClicked();

            void topsOriginalCHMClicked();
            void topsSmoothedCHMClicked();
            void topsTreetopsDatabaseClicked();

            void crownsCrownsDatabaseClicked();
            void crownsCrownsRasterClicked();
            void crownsTreetopsDatabaseClicked();
            void crownsSmoothedCHMClicked();

            void done();
        };

        class WorkerThread : public QThread {
        private:
            TreetopsForm *m_parent;

            void run() {
                
                using namespace geotools::treetops;
                using namespace geotools::treetops::config;
                
                try {
                    Treetops t;
                    t.setCallbacks(m_parent->m_callbacks);
                    const TreetopsConfig &config = m_parent->m_config;
                    const TreetopsCallbacks *cb = (TreetopsCallbacks *) m_parent->m_callbacks;
                    
                    int steps = ((int) config.doSmoothing) + ((int) config.doTops) + ((int) config.doCrowns);
                    int step = 0;
                    
                    if(cb)
                        cb->overallCallback(0.01);
                    
                    if(config.doSmoothing) {
                        t.smooth(config);
                        cb->overallCallback((float) ++step / steps);
                    }
                    
                    if(config.doTops) {
                        t.treetops(config);
                        cb->overallCallback((float) ++step / steps);
                    }
                    
                    if(config.doCrowns) {
                        t.treecrowns(config);
                        cb->overallCallback((float) ++step / steps);
                    }
                    
                    cb->overallCallback(1.0);
                    
                } catch (const std::exception &e) {
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


#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QMessageBox>
#include <QDir>
//#include <QHelpEngine>

#include "geotools.h"
#include "treetops.hpp"
#include "treetops_ui.hpp"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;
using namespace geotools::treetops;

TreetopsForm::TreetopsForm(QWidget *p) :
    QWidget(p),
    m_callbacks(nullptr) {
}

TreetopsForm::~TreetopsForm() {
    delete m_callbacks;
    if (m_workerThread) {
        m_workerThread->exit(0);
        delete m_workerThread;
    }
}

void TreetopsForm::setupUi(QWidget *form) {

    Ui::TreetopsForm::setupUi(form);

    g_loglevel(G_LOG_DEBUG);

    m_callbacks = new TreetopsCallbacks();

    m_workerThread = new WorkerThread();
    m_workerThread->init(this);
    
    m_form = form;
    m_last = QDir::home();

    connect(chkEnableSmoothing, SIGNAL(toggled(bool)), SLOT(doSmoothChanged(bool)));
    connect(chkEnableTops, SIGNAL(toggled(bool)), SLOT(doTopsChanged(bool)));
    connect(chkEnableCrowns, SIGNAL(toggled(bool)), SLOT(doCrownsChanged(bool)));

    connect(spnSmoothWindow, SIGNAL(valueChanged(int)), SLOT(smoothWindowSizeChanged(int)));
    connect(spnSmoothStdDev, SIGNAL(valueChanged(double)), SLOT(smoothStdDevChanged(double)));
    connect(txtSmoothOriginalCHM, SIGNAL(textChanged(QString)), SLOT(smoothOriginalCHMChanged(QString)));
    connect(txtSmoothSmoothedCHM, SIGNAL(textChanged(QString)), SLOT(smoothSmoothedCHMChanged(QString)));
    connect(btnSmoothOriginalCHM, SIGNAL(clicked()), SLOT(smoothOriginalCHMClicked()));
    connect(btnSmoothSmoothedCHM, SIGNAL(clicked()), SLOT(smoothSmoothedCHMClicked()));

    connect(spnTopsMinHeight, SIGNAL(valueChanged(double)), SLOT(topsMinHeightChanged(double)));
    connect(spnTopsWindowSize, SIGNAL(valueChanged(int)), SLOT(topsWindowSizeChanged(int)));
    connect(spnTopsTreetopsSRID, SIGNAL(valueChanged(int)), SLOT(topsTreetopsSRIDChanged(int)));
    connect(txtTopsOriginalCHM, SIGNAL(textChanged(QString)), SLOT(topsOriginalCHMChanged(QString)));
    connect(txtTopsSmoothedCHM, SIGNAL(textChanged(QString)), SLOT(topsSmoothedCHMChanged(QString)));
    connect(txtTopsTreetopsDatabase, SIGNAL(textChanged(QString)), SLOT(topsTreetopsDatabaseChanged(QString)));
    connect(btnTopsOriginalCHM, SIGNAL(clicked()), SLOT(topsOriginalCHMClicked()));
    connect(btnTopsSmoothedCHM, SIGNAL(clicked()), SLOT(topsSmoothedCHMClicked()));
    connect(btnTopsTreetopsDatabase, SIGNAL(clicked()), SLOT(topsTreetopsDatabaseClicked()));

    connect(spnCrownsRadius, SIGNAL(valueChanged(double)), SLOT(crownsRadiusChanged(double)));
    connect(spnCrownsHeightFraction, SIGNAL(valueChanged(double)), SLOT(crownsHeightFractionChanged(double)));
    connect(spnCrownsMinHeight, SIGNAL(valueChanged(double)), SLOT(crownsMinHeightChanged(double)));
    connect(txtCrownsSmoothedCHM, SIGNAL(textChanged(QString)), SLOT(crownsSmoothedCHMChanged(QString)));
    connect(txtCrownsTreetopsDatabase, SIGNAL(textChanged(QString)), SLOT(crownsTreetopsDatabaseChanged(QString)));
    connect(txtCrownsCrownsRaster, SIGNAL(textChanged(QString)), SLOT(crownsCrownsRasterChanged(QString)));
    connect(txtCrownsCrownsDatabase, SIGNAL(textChanged(QString)), SLOT(crownsCrownsDatabaseChanged(QString)));
    connect(btnCrownsSmoothedCHM, SIGNAL(clicked()), SLOT(crownsSmoothedCHMClicked()));
    connect(btnCrownsTreetopsDatabase, SIGNAL(clicked()), SLOT(crownsTreetopsDatabaseClicked()));
    connect(btnCrownsCrownsRaster, SIGNAL(clicked()), SLOT(crownsCrownsRasterClicked()));
    connect(btnCrownsCrownsDatabase, SIGNAL(clicked()), SLOT(crownsCrownsDatabaseClicked()));
    connect(btnTopsTreetopsSRID, SIGNAL(clicked()), SLOT(topsTreetopsSRIDClicked()));

    connect(btnExit, SIGNAL(clicked()), SLOT(exitClicked()));
    connect(btnRun, SIGNAL(clicked()), SLOT(runClicked()));
    connect(btnCancel, SIGNAL(clicked()), SLOT(cancelClicked()));
    connect(btnHelp, SIGNAL(clicked()), SLOT(helpClicked()));
    
    if (m_callbacks) {
        connect((TreetopsCallbacks *) m_callbacks, SIGNAL(stepProgress(int)), prgStep, SLOT(setValue(int)));
        connect((TreetopsCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
    }
    connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));

}

QString _str(const std::string &str) {
    return QString(str.c_str());
}

void TreetopsForm::topsTreetopsSRIDChanged(int srid) {
    m_config.srid = spnTopsTreetopsSRID->value();
    checkRun();
}

void TreetopsForm::topsTreetopsSRIDClicked() {
    CRSSelector cs(m_form);
    cs.enableVertical(false);
    cs.setHorizontalSRID(m_config.srid);
    if (cs.exec())
        spnTopsTreetopsSRID->setValue(cs.getHorizontalSRID());
}

void TreetopsForm::updateView() {
    /*
    txtSmoothInputFile->setText(_str(m_config.smoothInputFile));
    txtSmoothOutputFile->setText(_str(m_config.smoothOutputFile));
    txtTopsOriginalCHM->setText(_str(m_config.topsInputFile));
    txtTopsOutputFile->setText(_str(m_config.topsOutputFile));
    txtCrownsInputFile->setText(_str(m_config.crownsInputFile));
    txtCrownsTreetopsFile->setText(_str(m_config.crownsTreetopsFile));
    txtCrownsOutputRaster->setText(_str(m_config.crownsOutputRaster));
    txCrownsOutputVector->setText(_str(m_config.crownsOutputVector));
     */
}

std::string _getInputFile(QWidget *form, const std::string &title, QDir &path, const std::string &filter) {
    QString res = QFileDialog::getOpenFileName(form, _str(title), path.path(), _str(filter));
    path.setPath(res);
    return res.toStdString();
}

std::string _getOutputFile(QWidget *form, const std::string &title, QDir &path, const std::string &filter) {
    QString res = QFileDialog::getSaveFileName(form, _str(title), path.path(), _str(filter));
    path.setPath(res);
    return res.toStdString();
}

void TreetopsForm::smoothOriginalCHMClicked() {
    m_config.smoothOriginalCHM = _getInputFile(this, "CHM for Smoothing", m_last, "GeoTiff (*.tif *.tiff)");
    txtSmoothOriginalCHM->setText(_str(m_config.smoothOriginalCHM));
    checkRun();
    if (m_config.topsOriginalCHM.empty()) {
        m_config.topsOriginalCHM = m_config.smoothOriginalCHM;
        txtTopsOriginalCHM->setText(_str(m_config.topsOriginalCHM));
    }
    updateView();
}

void TreetopsForm::smoothSmoothedCHMClicked() {
    m_config.smoothSmoothedCHM = _getOutputFile(this, "Smoothed CHM", m_last, "GeoTiff (*.tif *.tiff)");
    txtSmoothSmoothedCHM->setText(_str(m_config.smoothSmoothedCHM));
    checkRun();
    if (m_config.topsSmoothedCHM.empty()) {
        m_config.topsSmoothedCHM = m_config.smoothSmoothedCHM;
        txtTopsSmoothedCHM->setText(_str(m_config.topsSmoothedCHM));
    }
    if (m_config.crownsSmoothedCHM.empty()) {
        m_config.crownsSmoothedCHM = m_config.smoothSmoothedCHM;
        txtCrownsSmoothedCHM->setText(_str(m_config.crownsSmoothedCHM));
    }
    updateView();
}

void TreetopsForm::topsSmoothedCHMClicked() {
    m_config.topsSmoothedCHM = _getInputFile(this, "Smoothed CHM for Treetops", m_last, "GeoTiff (*.tif *.tiff)");
    txtTopsSmoothedCHM->setText(_str(m_config.topsSmoothedCHM));
    checkRun();
    if (m_config.crownsSmoothedCHM.empty()) {
        m_config.crownsSmoothedCHM = m_config.topsSmoothedCHM;
        txtCrownsSmoothedCHM->setText(_str(m_config.crownsSmoothedCHM));
    }
    updateView();
}

void TreetopsForm::topsOriginalCHMClicked() {
    m_config.topsOriginalCHM = _getInputFile(this, "Original CHM for Treetops", m_last, "GeoTiff (*.tif *.tiff)");
    txtTopsOriginalCHM->setText(_str(m_config.topsOriginalCHM));
    checkRun();
    updateView();
}

void TreetopsForm::topsTreetopsDatabaseClicked() {
    m_config.topsTreetopsDatabase = _getOutputFile(this, "Treetops Database", m_last, "SQLite (*.sqlite)");
    txtTopsTreetopsDatabase->setText(_str(m_config.topsTreetopsDatabase));
    checkRun();
    if (m_config.crownsTreetopsDatabase.empty()) {
        m_config.crownsTreetopsDatabase = m_config.topsTreetopsDatabase;
        txtCrownsTreetopsDatabase->setText(_str(m_config.crownsTreetopsDatabase));
    }
    updateView();
}

void TreetopsForm::crownsSmoothedCHMClicked() {
    m_config.crownsSmoothedCHM = _getInputFile(this, "Smoothed CHM for Crown Delineation", m_last, "GeoTiff (*.tif *.tiff)");
    txtCrownsSmoothedCHM->setText(_str(m_config.crownsSmoothedCHM));
    checkRun();
}

void TreetopsForm::crownsTreetopsDatabaseClicked() {
    m_config.crownsTreetopsDatabase = _getInputFile(this, "Treetops Database", m_last, "SQLite (*.sqlite)");
    txtCrownsTreetopsDatabase->setText(_str(m_config.crownsTreetopsDatabase));
    checkRun();
}

void TreetopsForm::crownsCrownsRasterClicked() {
    m_config.crownsCrownsRaster = _getOutputFile(this, "Crowns Raster", m_last, "GeoTiff (*.tif *.tiff)");
    txtCrownsCrownsRaster->setText(_str(m_config.crownsCrownsRaster));
    checkRun();
}

void TreetopsForm::crownsCrownsDatabaseClicked() {
    m_config.crownsCrownsDatabase = _getOutputFile(this, "Crowns Database", m_last, "SQLite (*.sqlite)");
    txtCrownsCrownsDatabase->setText(_str(m_config.crownsCrownsDatabase));
    checkRun();
}

void TreetopsForm::doSmoothChanged(bool v) {
    m_config.doSmoothing = v;
    checkRun();
}

void TreetopsForm::doTopsChanged(bool v) {
    m_config.doTops = v;
    checkRun();
}

void TreetopsForm::doCrownsChanged(bool v) {
    m_config.doCrowns = v;
    checkRun();
}

void TreetopsForm::crownsRadiusChanged(double radius) {
    m_config.crownsRadius = radius;
    checkRun();
}

void TreetopsForm::crownsHeightFractionChanged(double frac) {
    m_config.crownsHeightFraction = frac;
    checkRun();
}

void TreetopsForm::crownsMinHeightChanged(double height) {
    m_config.crownsMinHeight = height;
    checkRun();
}

void TreetopsForm::crownsSmoothedCHMChanged(QString file) {
    m_config.crownsSmoothedCHM = file.toStdString();
    checkRun();
}

void TreetopsForm::crownsTreetopsDatabaseChanged(QString file) {
    m_config.crownsTreetopsDatabase = file.toStdString();
    checkRun();
}

void TreetopsForm::crownsCrownsRasterChanged(QString file) {
    m_config.crownsCrownsRaster = file.toStdString();
    checkRun();
}

void TreetopsForm::crownsCrownsDatabaseChanged(QString file) {
    m_config.crownsCrownsDatabase = file.toStdString();
    checkRun();
}

void TreetopsForm::topsMinHeightChanged(double height) {
    m_config.topsMinHeight = height;
    checkRun();
}

void TreetopsForm::topsWindowSizeChanged(int size) {
    m_config.topsWindowSize = size;
    checkRun();
}

void TreetopsForm::topsSmoothedCHMChanged(QString file) {
    m_config.topsSmoothedCHM = file.toStdString();
    checkRun();
}

void TreetopsForm::topsOriginalCHMChanged(QString file) {
    m_config.topsOriginalCHM = file.toStdString();
    checkRun();
}

void TreetopsForm::topsTreetopsDatabaseChanged(QString file) {
    m_config.topsTreetopsDatabase = file.toStdString();
    checkRun();
}

void TreetopsForm::smoothWindowSizeChanged(int size) {
    m_config.smoothWindowSize = size;
    checkRun();
}

void TreetopsForm::smoothStdDevChanged(double stdDev) {
    m_config.smoothStdDev = stdDev;
    checkRun();
}

void TreetopsForm::smoothOriginalCHMChanged(QString file) {
    m_config.smoothOriginalCHM = file.toStdString();
    checkRun();
}

void TreetopsForm::smoothSmoothedCHMChanged(QString file) {
    m_config.smoothSmoothedCHM = file.toStdString();
    checkRun();
}

/*
void TreetopsForm::crsConfigClicked() {
        CRSSelector cs(m_form);
        cs.setHorizontalSRID(m_hsrid);
        cs.setVerticalSRID(m_vsrid);
        if(cs.exec()) {
                m_vsrid = cs.getVerticalSRID();
                m_hsrid = cs.getHorizontalSRID();
                std::stringstream ss;
                ss << "epsg:" << m_hsrid;
                if(m_vsrid > 0)
                        ss << "+" << m_vsrid;
                txtCRSConfig->setText(QString(ss.str().c_str()));
        }
        checkRun();
}
 */

void TreetopsForm::runClicked() {

    if (m_workerThread->isRunning())
        return;

    btnRun->setEnabled(false);
    btnExit->setEnabled(false);
    btnCancel->setEnabled(true);

    m_workerThread->init(this);
    m_workerThread->start();

}

void TreetopsForm::done() {
    btnRun->setEnabled(true);
    btnExit->setEnabled(true);
    btnCancel->setEnabled(false);
    checkRun();
}

void TreetopsForm::exitClicked() {
    g_trace("quit");
    m_form->close();
}

void TreetopsForm::cancelClicked() {
    g_trace("cancel");
}

void TreetopsForm::helpClicked() {
    //QHelpEngine qh("help.qhc", this);
    g_trace("help");
}

void TreetopsForm::checkRun() {
    // TODO: Check runnable.
    btnRun->setEnabled(true);
}
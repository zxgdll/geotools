#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QMessageBox>
#include <QDir>

#include "geotools.h"
#include "treetops.hpp"
#include "treetops_ui.hpp"

using namespace geotools::ui;
using namespace geotools::trees;


TreetopsForm::TreetopsForm(QWidget *p) : 
	QWidget(p),
	m_callbacks(nullptr) {
}

TreetopsForm::~TreetopsForm() {
	delete m_callbacks;
	if(m_workerThread) {
		m_workerThread->exit(0);
		delete m_workerThread;
	}
}

void TreetopsForm::setupUi(QWidget *form) {

	Ui::TreetopsForm::setupUi(form);

	m_workerThread = new WorkerThread();
	m_workerThread->init(this);
	
	g_loglevel(G_LOG_DEBUG);

	m_form = form;
	m_last = QDir::home();

	connect(chkEnableSmoothing, SIGNAL(toggled(bool)), SLOT(doSmoothChanged(bool)));
	connect(chkEnableTops, SIGNAL(toggled(bool)), SLOT(doTopsChanged(bool)));
	connect(chkEnableCrowns, SIGNAL(toggled(bool)), SLOT(doCrownsChanged(bool)));

	connect(spnSmoothWindow, SIGNAL(valueChanged(int)), SLOT(smoothWindowSizeChanged(int)));
	connect(spnSmoothStdDev, SIGNAL(valueChanged(double)), SLOT(smoothStdDevChanged(double)));
	connect(txtSmoothInputFile, SIGNAL(textChanged(QString)), SLOT(smoothInputFileChanged(QString)));
	connect(txtSmoothOutputFile, SIGNAL(textChanged(QString)), SLOT(smoothOutputFileChanged(QString)));
	connect(btnSmoothInputFile, SIGNAL(clicked()), SLOT(smoothInputFileClicked()));
	connect(btnSmoothOutputFile, SIGNAL(clicked()), SLOT(smoothOutputFileClicked()));

	connect(spnTopsMinHeight, SIGNAL(valueChanged(double)), SLOT(topsMinHeightChanged(double)));
	connect(spnTopsWindowSize, SIGNAL(valueChanged(int)), SLOT(topsWindowSizeChanged(int)));
	connect(txtTopsInputFile, SIGNAL(textChanged(QString)), SLOT(topsInputFileChanged(QString)));
	connect(txtTopsOutputFile, SIGNAL(textChanged(QString)), SLOT(topsOutputFileChanged(QString)));
	connect(btnTopsInputFile, SIGNAL(clicked()), SLOT(topsInputFileClicked()));
	connect(btnTopsOutputFile, SIGNAL(clicked()), SLOT(topsOutputFileClicked()));

	connect(spnCrownsRadius, SIGNAL(valueChanged(double)), SLOT(crownsRadiusChanged(double)));
	connect(spnCrownsHeightFraction, SIGNAL(valueChanged(double)), SLOT(crownsHeightFractionChanged(double)));
	connect(spnCrownsMinHeight, SIGNAL(valueChanged(double)), SLOT(crownsMinHeightChanged(double)));
	connect(txtCrownsInputFile, SIGNAL(textChanged(QString)), SLOT(crownsInputFileChanged(QString)));
	connect(txtCrownsTreetopsFile, SIGNAL(textChanged(QString)), SLOT(crownsTreetopsFileChanged(QString)));
	connect(txtCrownsOutputRaster, SIGNAL(textChanged(QString)), SLOT(crownsOutputRasterChanged(QString)));
	connect(txtCrownsOutputVector, SIGNAL(textChanged(QString)), SLOT(crownsOutputVectorChanged(QString)));
	connect(btnCrownsInputFile, SIGNAL(clicked()), SLOT(crownsInputFileClicked()));
	connect(btnCrownsTreetopsFile, SIGNAL(clicked()), SLOT(crownsTreetopsFileClicked()));
	connect(btnCrownsOutputRaster, SIGNAL(clicked()), SLOT(crownsOutputRasterClicked()));
	connect(btnCrownsOutputVector, SIGNAL(clicked()), SLOT(crownsOutputVectorClicked()));

	connect(btnExit, SIGNAL(clicked()), SLOT(exitClicked()));
	connect(btnRun, SIGNAL(clicked()), SLOT(runClicked()));
	connect(btnCancel, SIGNAL(clicked()), SLOT(cancelClicked()));

	if(m_callbacks) {
		connect((TreetopsCallbacks *) m_callbacks, SIGNAL(fileProgress(int)), prgStep, SLOT(setValue(int)));
		connect((TreetopsCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
	}
	connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));

}

std::string _getInputFile(QWidget *form, const std::string &title, const QDir &path, const std::string &filter) {
	QString res = QFileDialog::getOpenFileName(form, QString(title.c_str()), path.path(), QString(filter.c_str()));
	return res.toStdString();
}

std::string _getOutputFile(QWidget *form, const std::string &title, const QDir &path, const std::string &filter) {
	QString res = QFileDialog::getSaveFileName(form, QString(title.c_str()), path.path(), QString(filter.c_str()));
	return res.toStdString();
}

void TreetopsForm::smoothInputFileClicked() {
	m_config.smoothInputFile = _getInputFile(this, "CHM for Smoothing", m_last, "GeoTiff (*.tif; *.tiff)");
	checkRun();
}

void TreetopsForm::smoothOutputFileClicked() {
	m_config.smoothOutputFile = _getOutputFile(this, "Smoothed CHM", m_last, "GeoTiff (*.tif; *.tiff)");
	checkRun();
}

void TreetopsForm::topsInputFileClicked() {
	m_config.smoothInputFile = _getInputFile(this, "Smoothed CHM for Treetops", m_last, "GeoTiff (*.tif; *.tiff)");
	checkRun();
}

void TreetopsForm::topsOutputFileClicked() {
	m_config.smoothOutputFile = _getOutputFile(this, "Treetops Database", m_last, "SQLite (*.sqlite)");
	checkRun();
}

void TreetopsForm::crownsInputFileClicked() {
	m_config.crownsInputFile = _getInputFile(this, "Smoothed CHM for Crown Delineation", m_last, "GeoTiff (*.tif; *.tiff)");
	checkRun();
}

void TreetopsForm::crownsTreetopsFileClicked() {
	m_config.crownsTreetopsFile = _getInputFile(this, "Treetops Database", m_last, "SQLite (*.sqlite)");
	checkRun();
}

void TreetopsForm::crownsOutputRasterClicked() {
	m_config.crownsOutputRaster = _getOutputFile(this, "Crowns Output Raster", m_last, "GeoTiff (*.tif; *.tiff)");
	checkRun();
}

void TreetopsForm::crownsOutputVectorClicked() {
	m_config.crownsOutputVector = _getOutputFile(this, "Crowns Output Database", m_last, "SQLite (*.sqlite)");
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

void TreetopsForm::crownsInputFileChanged(QString file) {
	m_config.crownsInputFile = file.toStdString();
	checkRun();
}

void TreetopsForm::crownsTreetopsFileChanged(QString file) {
	m_config.crownsTreetopsFile = file.toStdString();
	checkRun();
}

void TreetopsForm::crownsOutputRasterChanged(QString file) {
	m_config.crownsOutputRaster = file.toStdString();
	checkRun();
}

void TreetopsForm::crownsOutputVectorChanged(QString file) {
	m_config.crownsOutputVector = file.toStdString();
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

void TreetopsForm::topsInputFileChanged(QString file) {
	m_config.topsInputFile = file.toStdString();
	checkRun();
}

void TreetopsForm::topsOutputFileChanged(QString file) {
	m_config.topsOutputFile = file.toStdString();
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

void TreetopsForm::smoothInputFileChanged(QString file) {
	m_config.smoothInputFile = file.toStdString();
	checkRun();
}

void TreetopsForm::smoothOutputFileChanged(QString file) {
	m_config.smoothOutputFile = file.toStdString();
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

	if(m_workerThread->isRunning())
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

void TreetopsForm::checkRun() {
	// TODO: Check runnable.
	btnRun->setEnabled(true);
}
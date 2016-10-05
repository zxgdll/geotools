#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QMessageBox>

#include "geotools.h"
#include "lasgrid.hpp"
#include "lasgrid_ui.hpp"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;
using namespace geotools::las::lasgrid_config;

/**
 * Utility methods for status callbacks.
 */
 LasgridForm *callbackForm;

void fileCallback(float status) {
	callbackForm->setFileStatus(status);
}

void overallCallback(float status) {
	callbackForm->setOverallStatus(status);
}



LasgridForm::LasgridForm(QWidget *p) : 
	QWidget(p),
	m_vsrid(0), m_hsrid(0) {
}

LasgridForm::~LasgridForm() {
	m_last = nullptr;
}

void LasgridForm::setupUi(QWidget *form) {

	Ui::LasgridForm::setupUi(form);

	g_loglevel(G_LOG_DEBUG);

	m_form = form;
	m_last = new QDir(QDir::home());

	m_radius = defaultRadius;
	m_resolution = defaultResolution;
	m_angleLimit = defaultAngleLimit;
	spnResolution->setValue(m_resolution);
	spnRadius->setValue(m_radius);
	spnMaxAngle->setValue(m_angleLimit);

	m_type = defaultType;
	int i = 0;
	int defaultIdx = -1;
	for(auto it = types.begin(); it != types.end(); ++it) {
		cboType->addItem(QString::fromStdString(it->first), QVariant(it->second));
		if(it->second == m_type)
			defaultIdx = i;
		++i;
	}
	cboType->setCurrentIndex(defaultIdx);

	m_attribute = defaultAttribute;
	i = 0;
	defaultIdx = -1;
	for(auto it = attributes.begin(); it != attributes.end(); ++it) {
		cboAttribute->addItem(QString::fromStdString(it->first), QVariant(it->second));
		if(it->second == m_attribute)
			defaultIdx = i;
		++i;
	}
	cboAttribute->setCurrentIndex(defaultIdx);

	for(i = 0; i < 255; ++i) {
		QString str;
		str.setNum(i);
		QListWidgetItem *item = new QListWidgetItem(str, lstClasses);
		item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
		item->setCheckState(defaultClasses.count(i) > 0 ? Qt::Checked : Qt::Unchecked);
		if(defaultClasses.count(i) > 0)
			m_classes.insert(i);
		lstClasses->addItem(item);
	}

	connect(btnSelectFiles, SIGNAL(clicked()), this, SLOT(selectFilesClicked()));
	connect(btnRemoveSelected, SIGNAL(clicked()), this, SLOT(removeFilesClicked()));
	connect(btnClearFiles, SIGNAL(clicked()), this, SLOT(clearFilesClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
	connect(btnDestFile, SIGNAL(clicked()), this, SLOT(destFileClicked()));
	connect(lstFiles, SIGNAL(itemSelectionChanged()), this, SLOT(fileListSelectionChanged()));
	connect(btnCRSConfig, SIGNAL(clicked()), this, SLOT(crsConfigClicked()));
	connect(cboType, SIGNAL(currentIndexChanged(int)), this, SLOT(typeSelected(int)));
	connect(spnResolution, SIGNAL(valueChanged(double)), this, SLOT(resolutionChanged(double)));
	connect(chkSnapToGrid, SIGNAL(toggled(bool)), this, SLOT(snapToGridChanged(bool)));
	connect(spnRadius, SIGNAL(valueChanged(double)), this, SLOT(radiusChanged(double)));
	connect(cboAttribute, SIGNAL(currentIndexChanged(int)), this, SLOT(attributeSelected(int)));
	connect(spnQuantile, SIGNAL(valueChanged(int)), SLOT(quantileChanged(int)));
	connect(spnQuantiles, SIGNAL(valueChanged(int)), SLOT(quantilesChanged(int)));
	connect(spnMaxAngle, SIGNAL(valueChanged(int)), SLOT(maxAngleChanged(int)));
	connect(this, SIGNAL(fileProgress(int)), prgFile, SLOT(setValue(int)));
	connect(this, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));

}

void LasgridForm::maxAngleChanged(int q) {
	m_angleLimit = q;
	checkRun();
}

void LasgridForm::quantileChanged(int q) {
	m_quantile = q;
	checkRun();
}

void LasgridForm::quantilesChanged(int q) {
	m_quantiles = q;
	checkRun();
}

void LasgridForm::attributeSelected(int index) {
	std::string att = cboAttribute->itemText(index).toStdString();
	m_attribute = attributes[att];
	checkRun();
}

void LasgridForm::radiusChanged(double radius) {
	m_radius = radius;
	checkRun();
}

void LasgridForm::snapToGridChanged(bool state) {
	m_snap = state;
	checkRun();
}

void LasgridForm::resolutionChanged(double res) {
	m_resolution = res;
	checkRun();
}

void LasgridForm::typeSelected(int index) {
	std::string type = cboType->itemText(index).toStdString();
	m_type = types[type];
	bool on = m_type == TYPE_QUANTILE;
	spnQuantile->setEnabled(on);
	spnQuantiles->setEnabled(on);
	checkRun();
}

void LasgridForm::crsConfigClicked() {
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

void LasgridForm::fileListSelectionChanged() {
	updateFileButtons();
	checkRun();
}

void LasgridForm::destFileClicked() {
	QFileDialog d(m_form);
	d.setDirectory(*m_last);
	//d.setFileMode(QFileDialog::ExistingFile);
	d.setNameFilter(QString("GeoTiff Files (*.tif)"));
	if(d.exec()) {
		m_destFile = d.selectedFiles()[0].toStdString();
	} else {
		m_destFile = "";
	}
	txtDestFile->setText(QString(m_destFile.c_str()));
	checkRun();
}

void LasgridForm::setFileStatus(float status) {
	emit fileProgress((int) std::round(status * 100));
}

void LasgridForm::setOverallStatus(float status) {
	emit overallProgress((int) std::round(status * 100));
}

void LasgridForm::runClicked() {
	g_trace("run");

	using namespace geotools::las;
	using namespace geotools::util;

	Bounds bounds;

	m_angleLimit = 100;

	// TODO: Bounds
	// TODO: Angle limit
	// TODO: Compound CRS
	try {
		btnRun->setEnabled(false);
		btnCancel->setEnabled(false);
		callbackForm = this;
		LasGrid lg;
		lg.setFileCallback(fileCallback);
		lg.setOverallCallback(overallCallback);
		lg.lasgrid(m_destFile, m_lasFiles, m_classes, m_hsrid, m_attribute, m_type, m_radius, m_resolution, 
			bounds, m_angleLimit, m_fill, m_snap);
		checkRun();
		btnCancel->setEnabled(true);
	} catch(const std::exception &e) {
		QMessageBox err(this);
		err.setText("Error");
		err.setInformativeText(QString(e.what()));
		err.exec();
		checkRun();
		btnCancel->setEnabled(true);
	}
}

void LasgridForm::cancelClicked() {
	g_trace("quit");
	m_form->close();
}

void LasgridForm::updateFileList() {
	while(lstFiles->count())
		lstFiles->takeItem(0);
	for(unsigned int i = 0; i < m_lasFiles.size(); ++i)
		lstFiles->addItem(QString(m_lasFiles[i].c_str()));
	updateFileButtons();
	checkRun();
}

void LasgridForm::updateFileButtons() {
	btnClearFiles->setEnabled(lstFiles->count() > 0);
	btnRemoveSelected->setEnabled(lstFiles->selectedItems().size() > 0);
}

void LasgridForm::removeFilesClicked() {
	std::vector<std::string> lst;
	for(int i = 0; i < lstFiles->count(); ++i) {
		QListWidgetItem *item = lstFiles->item(i);
		if(!item->isSelected())
			lst.push_back(m_lasFiles[i]);
	}
	m_lasFiles.clear();
	m_lasFiles.assign(lst.begin(), lst.end());
	updateFileList();
	checkRun();
}

void LasgridForm::clearFilesClicked() {
	m_lasFiles.clear();
	updateFileList();
	checkRun();
}

void LasgridForm::selectFilesClicked() {
	QFileDialog d(m_form);
	d.setDirectory(*m_last);
	d.setFileMode(QFileDialog::ExistingFiles);
	d.setNameFilter(QString("LAS Files (*.las)"));
	if(d.exec()) {
		QStringList files = d.selectedFiles();
		*m_last = d.directory();
		std::set<std::string> tmp(m_lasFiles.begin(), m_lasFiles.end());
		for(int i = 0; i < files.size(); ++i)
			tmp.insert(files[i].toStdString());
		m_lasFiles.clear();
		m_lasFiles.assign(tmp.begin(), tmp.end());
	}
	updateFileList();
	checkRun();
}

void LasgridForm::checkRun() {
	// TODO: Check runnable.
	btnRun->setEnabled(true);
}
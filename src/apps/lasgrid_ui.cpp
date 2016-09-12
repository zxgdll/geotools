#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>

#include "geotools.h"
#include "lasgrid.hpp"
#include "lasgrid_ui.hpp"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;
using namespace geotools::las::lasgrid_config;

LasgridForm::LasgridForm(QWidget *p) : QWidget(p) {
}

LasgridForm::~LasgridForm() {
	m_last = nullptr;
}

void LasgridForm::setupUi(QWidget *form) {

	Ui::LasgridForm::setupUi(form);

	g_loglevel(G_LOG_TRACE);

	m_form = form;
	m_last = new QDir(QDir::home());

	spnResolution->setValue(defaultResolution);
	spnRadius->setValue(defaultRadius);

	int i = 0;
	int defaultIdx = -1;
	for(auto it = types.begin(); it != types.end(); ++it) {
		cboType->addItem(QString::fromStdString(it->first), QVariant(it->second));
		if(it->second == defaultType)
			defaultIdx = i;
		++i;
	}
	cboType->setCurrentIndex(defaultIdx);

	i = 0;
	defaultIdx = -1;
	for(auto it = attributes.begin(); it != attributes.end(); ++it) {
		cboAttribute->addItem(QString::fromStdString(it->first), QVariant(it->second));
		if(it->second == defaultAttribute)
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
}

void LasgridForm::crsConfigClicked() {
	CRSSelector cs(m_form);
	if(cs.exec());
}

void LasgridForm::fileListSelectionChanged() {
	updateFileButtons();
}

void LasgridForm::destFileClicked() {
	QFileDialog d(m_form);
	d.setDirectory(*m_last);
	d.setFileMode(QFileDialog::ExistingFile);
	d.setNameFilter(QString("GeoTiff Files (*.tif)"));
	if(d.exec()) {
		m_destFile = d.selectedFiles()[0].toStdString();
	} else {
		m_destFile = "";
	}
	txtDestFile->setText(QString(m_destFile.c_str()));
}

void LasgridForm::runClicked() {
	g_trace("run");

	using namespace geotools::las;
	using namespace geotools::util;

	Bounds bounds;

	lasgrid(m_destFile, m_lasFiles, m_classes, m_crs, m_attribute, m_type, m_radius, m_resolution, 
		bounds, m_angleLimit, m_fill, m_snap);
}

void LasgridForm::cancelClicked() {
	g_trace("quit");
	m_form->close();
}

void LasgridForm::updateFileList() {
	while(lstFiles->count())
		lstFiles->takeItem(0);
	for(int i = 0; i < m_lasFiles.size(); ++i)
		lstFiles->addItem(QString(m_lasFiles[i].c_str()));
	updateFileButtons();
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
}

void LasgridForm::clearFilesClicked() {
	m_lasFiles.clear();
	updateFileList();
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
}

#include "moc_lasgrid_ui.cpp"

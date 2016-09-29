#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>

#include "geotools.h"
#include "mosaic.hpp"
#include "mosaic_ui.hpp"

using namespace geotools::ui;

/**
 * Utility methods for status callbacks.
 */
MosaicForm *form;

void fileCallback(float status) {
	form->setFileStatus(status);
}

void overallCallback(float status) {
	form->setOverallStatus(status);
}



MosaicForm::MosaicForm(QWidget *p) : 
	QWidget(p),
	m_tileSize(2048),
	m_distance(100.0),
	m_threads(1) {
}

MosaicForm::~MosaicForm() {
	m_last = nullptr;
}

void MosaicForm::setupUi(QWidget *form) {

	Ui::MosaicForm::setupUi(form);

	g_loglevel(G_LOG_TRACE);

	m_form = form;
	m_last = new QDir(QDir::home());

	// TODO: Set from defaults -- see lasgrid_config
	spnThreads->setValue(m_threads);
	spnDistance->setValue(m_distance);
	spnTileSize->setValue(m_tileSize);

	connect(btnSelectFiles, SIGNAL(clicked()), this, SLOT(selectFilesClicked()));
	connect(btnRemoveSelected, SIGNAL(clicked()), this, SLOT(removeFilesClicked()));
	connect(btnClearFiles, SIGNAL(clicked()), this, SLOT(clearFilesClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
	connect(btnDestFile, SIGNAL(clicked()), this, SLOT(destFileClicked()));
	connect(lstFiles, SIGNAL(itemSelectionChanged()), this, SLOT(fileListSelectionChanged()));
	connect(spnDistance, SIGNAL(valueChanged(double)), this, SLOT(distanceChanged(double)));
	connect(spnTileSize, SIGNAL(valueChanged(int)), this, SLOT(tileSizeChanged(int)));
	connect(spnThreads, SIGNAL(valueChanged(int)), this, SLOT(threadsChanged(int)));
	connect(btnUp, SIGNAL(clicked()), this, SLOT(upClicked()));
	connect(btnDown, SIGNAL(clicked()), this, SLOT(downClicked()));
}

void MosaicForm::threadsChanged(int threads) {
	m_threads = threads;
	checkRun();
}

void MosaicForm::threadsChanged(int tileSize) {
	m_tileSize = tileSize;
	checkRun();
}

void MosaicForm::distanceChanged(double distance) {
	m_distance = distance;
	checkRun();
}

void MosaicForm::fileListSelectionChanged() {
	updateFileButtons();
	checkRun();
}

void MosaicForm::destFileClicked() {
	QFileDialog d(m_form);
	d.setDirectory(*m_last);
	d.setNameFilter(QString("GeoTiff Files (*.tif)"));
	if(d.exec()) {
		m_destFile = d.selectedFiles()[0].toStdString();
	} else {
		m_destFile = "";
	}
	txtDestFile->setText(QString(m_destFile.c_str()));
	checkRun();
}

void MosaicForm::setFileStatus(float status) {
	prgFile->setValue((int) status);
}

void MosaicForm::setOverallStatus(float status) {
	prgOverall->setValue((int) status);
}

void MosaicForm::runClicked() {
	g_trace("run");

//	try {

//	} catch(const std::exception &e) {
//		std::cerr << e.what() << std::endl;
//	}
}

void MosaicForm::cancelClicked() {
	m_form->close();
}

void MosaicForm::updateFileList() {
	while(lstFiles->count())
		lstFiles->takeItem(0);
	for(int i = 0; i < m_lasFiles.size(); ++i)
		lstFiles->addItem(QString(m_lasFiles[i].c_str()));
	updateFileButtons();
	checkRun();
}

void MosaicForm::updateFileButtons() {
	bool selected = lstFiles->selectedItems().size() > 0;
	btnClearFiles->setEnabled(lstFiles->count() > 0);
	btnRemoveSelected->setEnabled(selected);
	btnUp->setEnabled(selected && lstFiles->currentRow() > 0);
	btnDown->setEnabled(selected && lstFiles->currentRow() < lstFiles->count() - 1);
}

void MosaicForm::removeFilesClicked() {
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

void MosaicForm::clearFilesClicked() {
	m_lasFiles.clear();
	updateFileList();
	checkRun();
}

void MosaicForm::selectFilesClicked() {
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

void MosaicForm::checkRun() {
	// TODO: Check runnable.
	btnRun->setEnabled(true);
}
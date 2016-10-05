#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>

#include "omp.h"

#include "geotools.h"
#include "mosaic.hpp"
#include "mosaic_ui.hpp"

using namespace geotools::ui;
using namespace geotools::raster;


MosaicForm::MosaicForm(QWidget *p) : 
	QWidget(p),
	m_tileSize(2048),
	m_distance(100.0),
	m_threads(1),
	m_callbacks(new MosaicCallbacks()) {
}

MosaicForm::~MosaicForm() {
	m_last = nullptr;
	delete m_callbacks;
}

void MosaicForm::setupUi(QWidget *form) {

	Ui::MosaicForm::setupUi(form);

	g_loglevel(G_LOG_TRACE);

	m_form = form;
	m_last = new QDir(QDir::home());

	// TODO: Set from defaults -- see lasgrid_config
	spnThreads->setValue(m_threads);
	spnThreads->setMaximum(g_max(1, omp_get_num_procs()));

	QString q;
	q.setNum(m_distance);
	txtDistance->setText(q);
	q.setNum(m_tileSize);
	txtTileSize->setText(q);

	connect(btnSelectFiles, SIGNAL(clicked()), this, SLOT(selectFilesClicked()));
	connect(btnRemoveSelected, SIGNAL(clicked()), this, SLOT(removeFilesClicked()));
	connect(btnClearFiles, SIGNAL(clicked()), this, SLOT(clearFilesClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
	connect(btnDestFile, SIGNAL(clicked()), this, SLOT(destFileClicked()));
	connect(lstFiles, SIGNAL(itemSelectionChanged()), this, SLOT(fileListSelectionChanged()));
	connect(txtDistance, SIGNAL(textChanged(QString)), this, SLOT(distanceChanged(QString)));
	connect(txtTileSize, SIGNAL(textChanged(QString)), this, SLOT(tileSizeChanged(QString)));
	connect(spnThreads, SIGNAL(valueChanged(int)), this, SLOT(threadsChanged(int)));
	connect(btnOrderUp, SIGNAL(clicked()), this, SLOT(upClicked()));
	connect(btnOrderDown, SIGNAL(clicked()), this, SLOT(downClicked()));
	connect((MosaicCallbacks *) m_callbacks, SIGNAL(fileProgress(int)), prgFile, SLOT(setValue(int)));
	connect((MosaicCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
	connect(&m_workerThread, SIGNAL(finished()), this, SLOT(done()));

}

void MosaicForm::upClicked() {
	int row = lstFiles->currentRow();
	std::string f = m_tifFiles[row];
	m_tifFiles[row] = m_tifFiles[row - 1];
	m_tifFiles[row - 1] = f;
	updateFileList();
}

void MosaicForm::downClicked() {
	int row = lstFiles->currentRow();
	std::string f = m_tifFiles[row];
	m_tifFiles[row] = m_tifFiles[row + 1];
	m_tifFiles[row + 1] = f;
	updateFileList();
}

void MosaicForm::threadsChanged(int threads) {
	m_threads = threads;
	checkRun();
}

void MosaicForm::tileSizeChanged(QString data) {
	m_tileSize = data.toInt();
	checkRun();
}

void MosaicForm::distanceChanged(QString data) {
	m_distance = data.toDouble();
	checkRun();
}

void MosaicForm::fileListSelectionChanged() {
	updateFileButtons();
	checkRun();
}

void MosaicForm::destFileClicked() {
	QString res = QFileDialog::getSaveFileName(this, "Save File", m_last->path(), "GeoTiff (*.tif *.tiff)");
	m_destFile = res.toStdString();
	m_last->setPath(res);
	txtDestFile->setText(QString(m_destFile.c_str()));
	checkRun();
}

void MosaicForm::runClicked() {
	if(m_workerThread.isRunning())
		return;
	btnRun->setEnabled(false);
	btnCancel->setEnabled(false);
	m_workerThread.init(this, m_callbacks, m_tifFiles, m_destFile, m_distance, m_tileSize, m_threads);
	m_workerThread.start();
}

void MosaicForm::done() {
	checkRun();
	btnCancel->setEnabled(true);
}

void MosaicForm::cancelClicked() {
	m_form->close();
}

void MosaicForm::updateFileList() {
	while(lstFiles->count())
		lstFiles->takeItem(0);
	for(unsigned int i = 0; i < m_tifFiles.size(); ++i)
		lstFiles->addItem(QString(m_tifFiles[i].c_str()));
	updateFileButtons();
	checkRun();
}

void MosaicForm::updateFileButtons() {
	bool selected = lstFiles->selectedItems().size() > 0;
	btnClearFiles->setEnabled(lstFiles->count() > 0);
	btnRemoveSelected->setEnabled(selected);
	btnOrderUp->setEnabled(selected && lstFiles->currentRow() > 0);
	btnOrderDown->setEnabled(selected && lstFiles->currentRow() < lstFiles->count() - 1);
}

void MosaicForm::removeFilesClicked() {
	std::vector<std::string> lst;
	for(int i = 0; i < lstFiles->count(); ++i) {
		QListWidgetItem *item = lstFiles->item(i);
		if(!item->isSelected())
			lst.push_back(m_tifFiles[i]);
	}
	m_tifFiles.clear();
	m_tifFiles.assign(lst.begin(), lst.end());
	updateFileList();
	checkRun();
}

void MosaicForm::clearFilesClicked() {
	m_tifFiles.clear();
	updateFileList();
	checkRun();
}

void MosaicForm::selectFilesClicked() {
	QFileDialog d(m_form);
	d.setDirectory(*m_last);
	d.setFileMode(QFileDialog::ExistingFiles);
	d.setNameFilter(QString("TIFF Files (*.tif)"));
	if(d.exec()) {
		QStringList files = d.selectedFiles();
		*m_last = d.directory();
		std::set<std::string> tmp(m_tifFiles.begin(), m_tifFiles.end());
		for(int i = 0; i < files.size(); ++i)
			tmp.insert(files[i].toStdString());
		m_tifFiles.clear();
		m_tifFiles.assign(tmp.begin(), tmp.end());
	}
	updateFileList();
	checkRun();
}

void MosaicForm::checkRun() {
	bool enabled = true;
	if(m_threads <= 0)
		enabled = false;
	if(m_distance <= 0.0)
		enabled = false;
	if(m_tifFiles.size() < 2)
		enabled = false;
	btnRun->setEnabled(enabled);
}
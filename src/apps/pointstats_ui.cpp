#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QMessageBox>
#include <QSettings>

#include "geotools.h"
#include "pointstats.hpp"
#include "pointstats_ui.hpp"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;
using namespace geotools::point;
using namespace geotools::point::pointstats_config;

QSettings _settings("PointStats", "Geotools");
QString _last_dir("last_dir");

void WorkerThread::run() {
	try {
		m_error = "";
		geotools::point::PointStats l;
		
		geotools::point::PointStatsConfig config;
		config.dstFile = m_parent->m_destFile;
		config.sourceFiles = m_parent->m_sourceFiles;
		config.classes = m_parent->m_classes;
		config.hsrid = m_parent->m_hsrid;
		config.attribute = m_parent->m_attribute;
		config.type = m_parent->m_type;
		config.resolution = m_parent->m_resolution;
		config.bounds = m_bounds;
		config.angleLimit = m_parent->m_angleLimit;
		config.fill = m_parent->m_fill;
		config.snap = m_parent->m_snap;
		config.threads = m_parent->m_threads;
		config.gapFractionType = m_parent->m_gapFunction;

		l.pointstats(config, m_parent->m_callbacks);
	} catch(const std::exception &e) {
		m_error = e.what();
	}
}

bool WorkerThread::hasError() {
	return !m_error.empty();
}

std::string WorkerThread::getError() {
	return m_error;
}

PointStatsForm::PointStatsForm(QWidget *p) : 
	QWidget(p),
	m_vsrid(0), m_hsrid(0) {
}

PointStatsForm::~PointStatsForm() {
	delete m_callbacks;
	if(m_workerThread) {
		m_workerThread->exit(0);
		delete m_workerThread;
	}
}

void PointStatsForm::setupUi(QWidget *form) {

	Ui::PointStatsForm::setupUi(form);

	m_form = form;
	if(_settings.contains(_last_dir)) {
		m_last.setPath(_settings.value(_last_dir).toString());
	} else {
		m_last = QDir::home();
	}

	m_resolution = defaultResolution;
	m_angleLimit = defaultAngleLimit;
	spnResolution->setValue(m_resolution);
	spnMaxAngle->setValue(m_angleLimit);

	m_workerThread = new WorkerThread();
	m_callbacks = new PointStatsCallbacks();

	spnThreads->setValue(m_threads);
	spnThreads->setMaximum(g_max(1, omp_get_num_procs()));

	m_type = defaultType;
	int i = 0;
	int defaultIdx = -1;
	for(const auto &it : types) {
		cboType->addItem(QString::fromStdString(it.first), QVariant(it.second));
		if(it.second == m_type)
			defaultIdx = i;
		++i;
	}
	cboType->setCurrentIndex(defaultIdx);

	m_gapFunction = defaultGapFraction;
	i = 0;
	defaultIdx = -1;
	for(const auto &it : gapFractionTypes) {
		cboGapFunction->addItem(QString::fromStdString(it.first), QVariant(it.second));
		if(it.second == m_gapFunction)
			defaultIdx = i;
		++i;
	}
	cboGapFunction->setCurrentIndex(defaultIdx);

	m_attribute = defaultAttribute;
	i = 0;
	defaultIdx = -1;
	for(const auto &it : attributes) {
		cboAttribute->addItem(QString::fromStdString(it.first), QVariant(it.second));
		if(it.second == m_attribute)
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

	chkSnapToGrid->setCheckState(defaultSnapToGrid ? Qt::Checked : Qt::Unchecked);

	connect(btnSelectFiles, SIGNAL(clicked()), this, SLOT(selectFilesClicked()));
	connect(btnRemoveSelected, SIGNAL(clicked()), this, SLOT(removeFilesClicked()));
	connect(btnClearFiles, SIGNAL(clicked()), this, SLOT(clearFilesClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(runClicked()));
	connect(btnDestFile, SIGNAL(clicked()), this, SLOT(destFileClicked()));
	connect(lstFiles, SIGNAL(itemSelectionChanged()), this, SLOT(fileListSelectionChanged()));
	connect(lstClasses, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(classItemClicked(QListWidgetItem*)));
	connect(btnCRSConfig, SIGNAL(clicked()), this, SLOT(crsConfigClicked()));
	connect(cboType, SIGNAL(currentIndexChanged(int)), this, SLOT(typeSelected(int)));
	connect(spnResolution, SIGNAL(valueChanged(double)), this, SLOT(resolutionChanged(double)));
	connect(chkSnapToGrid, SIGNAL(toggled(bool)), this, SLOT(snapToGridChanged(bool)));
	connect(cboAttribute, SIGNAL(currentIndexChanged(int)), this, SLOT(attributeSelected(int)));
	connect(cboGapFunction, SIGNAL(currentIndexChanged(int)), this, SLOT(gapFunctionSelected(int)));
	connect(spnThreads, SIGNAL(valueChanged(int)), SLOT(threadsChanged(int)));
	connect(spnQuantile, SIGNAL(valueChanged(int)), SLOT(quantileChanged(int)));
	connect(spnQuantiles, SIGNAL(valueChanged(int)), SLOT(quantilesChanged(int)));
	connect(spnMaxAngle, SIGNAL(valueChanged(int)), SLOT(maxAngleChanged(int)));
	connect((PointStatsCallbacks *) m_callbacks, SIGNAL(overallProgress(int)), prgOverall, SLOT(setValue(int)));
	connect(m_workerThread, SIGNAL(finished()), this, SLOT(done()));

	updateTypeUi();
}

void PointStatsForm::classItemClicked(QListWidgetItem *item) {
	unsigned char c = (unsigned char) item->text().toUShort();
	if(item->checkState() == Qt::Checked) {
		m_classes.insert(c);
	} else {
		m_classes.erase(c);
	}
}

void PointStatsForm::threadsChanged(int threads) {
	m_threads = threads;
	checkRun();
}

void PointStatsForm::maxAngleChanged(int q) {
	m_angleLimit = q;
	checkRun();
}

void PointStatsForm::quantileChanged(int q) {
	m_quantile = q;
	checkRun();
}

void PointStatsForm::quantilesChanged(int q) {
	m_quantiles = q;
	checkRun();
}

void PointStatsForm::attributeSelected(int index) {
	std::string att = cboAttribute->itemText(index).toStdString();
	m_attribute = attributes[att];
	checkRun();
}

void PointStatsForm::gapFunctionSelected(int index) {
	std::string gap = cboGapFunction->itemText(index).toStdString();
	m_gapFunction = gapFractionTypes[gap];
	checkRun();
}

void PointStatsForm::snapToGridChanged(bool state) {
	m_snap = state;
	checkRun();
}

void PointStatsForm::resolutionChanged(double res) {
	m_resolution = res;
	checkRun();
}

void PointStatsForm::updateTypeUi() {
	// TODO: See state machine framework
	spnQuantile->setVisible(m_type == TYPE_QUANTILE);
	spnQuantiles->setVisible(m_type == TYPE_QUANTILE);
	lblQuantile->setVisible(m_type == TYPE_QUANTILE);
	lblQuantiles->setVisible(m_type == TYPE_QUANTILE);
	cboGapFunction->setVisible(m_type == TYPE_GAP_FRACTION);
	lblGapFunction->setVisible(m_type == TYPE_GAP_FRACTION);
	lblAttribute->setVisible(m_type != TYPE_GAP_FRACTION);
	cboAttribute->setVisible(m_type != TYPE_GAP_FRACTION);
}

void PointStatsForm::typeSelected(int index) {
	std::string type = cboType->itemText(index).toStdString();
	m_type = types[type];
	updateTypeUi();
	checkRun();
}

void PointStatsForm::crsConfigClicked() {
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

void PointStatsForm::fileListSelectionChanged() {
	updateFileButtons();
	checkRun();
}

void PointStatsForm::destFileClicked() {
	QString res = QFileDialog::getSaveFileName(this, "Save File", m_last.path(), "GeoTiff (*.tif *.tiff)");
	m_destFile = res.toStdString();
	m_last.setPath(res);
	_settings.setValue(_last_dir, m_last.path());
	txtDestFile->setText(res);
	checkRun();
}

void PointStatsForm::runClicked() {

	if(m_workerThread->isRunning())
		return;

	btnRun->setEnabled(false);
	btnCancel->setEnabled(false);

	// TODO: Bounds
	Bounds bounds;
	m_workerThread->init(this, bounds);
	m_workerThread->start();

}

void PointStatsForm::done() {
	if(m_workerThread->hasError()) {
		QMessageBox err((QWidget *) this);
		err.setText("Error");
		err.setInformativeText(QString(m_workerThread->getError().c_str()));
		err.exec();
	}
	checkRun();
	btnCancel->setEnabled(true);
}

void PointStatsForm::cancelClicked() {
	g_trace("quit");
	m_form->close();
}

void PointStatsForm::updateFileList() {
	while(lstFiles->count())
		lstFiles->takeItem(0);
	for(const std::string &file : m_sourceFiles)
		lstFiles->addItem(QString(file.c_str()));
	updateFileButtons();
	checkRun();
}

void PointStatsForm::updateFileButtons() {
	btnClearFiles->setEnabled(lstFiles->count() > 0);
	btnRemoveSelected->setEnabled(lstFiles->selectedItems().size() > 0);
}

void PointStatsForm::removeFilesClicked() {
	std::vector<std::string> lst;
	unsigned int i = 0;
	for(const std::string &file : m_sourceFiles) {
		QListWidgetItem *item = lstFiles->item(i);
		if(!item->isSelected())
			lst.push_back(file);
		++i;
	}
	m_sourceFiles.clear();
	m_sourceFiles.assign(lst.begin(), lst.end());
	updateFileList();
	checkRun();
}

void PointStatsForm::clearFilesClicked() {
	m_sourceFiles.clear();
	updateFileList();
	checkRun();
}

void PointStatsForm::selectFilesClicked() {
	QFileDialog d(m_form);
	d.setDirectory(m_last);
	d.setFileMode(QFileDialog::ExistingFiles);
	d.setNameFilter(QString("LAS Files (*.las)"));
	if(d.exec()) {
		QStringList files = d.selectedFiles();
		m_last = d.directory();
		_settings.setValue(_last_dir, m_last.path());
		std::set<std::string> tmp(m_sourceFiles.begin(), m_sourceFiles.end());
		for(int i = 0; i < files.size(); ++i)
			tmp.insert(files[i].toStdString());
		m_sourceFiles.clear();
		m_sourceFiles.assign(tmp.begin(), tmp.end());
	}
	updateFileList();
	checkRun();
}

void PointStatsForm::checkRun() {
	// TODO: Check runnable.
	btnRun->setEnabled(true);
}
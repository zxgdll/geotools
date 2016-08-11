#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>

#include "lasgrid.hpp"
#include "lasgrid_ui.hpp"

using namespace geotools::ui;
using namespace geotools::las::lasgrid_config;

LasgridForm::LasgridForm(QWidget *p) : QObject(p) {
}

void LasgridForm::setupUi(QWidget *form) {

	Ui::LasgridForm::setupUi(form);

	m_form = form;

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

	connect(btnSelectFiles, SIGNAL(clicked()), this, SLOT(selectFileClicked()));
}

void LasgridForm::selectFileClicked() {
	QString file = QFileDialog::getOpenFileName(m_form, QString("LAS File"), QString("/"), QString("LAS Files (*.las)"));
}

#include "moc_lasgrid_ui.cpp"

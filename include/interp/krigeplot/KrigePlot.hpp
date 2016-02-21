/********************************************************************************
** Form generated from reading UI file 'variogramyJ7917.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef VARIOGRAMYJ7917_H
#define VARIOGRAMYJ7917_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QHeaderView>
#include "qwt_plot.h"
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>

#include <list>
#include "interp/InterpPoint.hpp"

QT_BEGIN_NAMESPACE

class Ui_dlgSemivariogram
{
public:
    QDialogButtonBox *buttonBox;
    QwtPlot *qwtVariogram;

    void setupUi(QDialog *dlgSemivariogram)
    {
        if (dlgSemivariogram->objectName().isEmpty())
            dlgSemivariogram->setObjectName(QString::fromUtf8("dlgSemivariogram"));
        dlgSemivariogram->resize(736, 439);
        buttonBox = new QDialogButtonBox(dlgSemivariogram);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(380, 400, 341, 32));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        qwtVariogram = new QwtPlot(dlgSemivariogram);
        qwtVariogram->setObjectName(QString::fromUtf8("qwtVariogram"));
        qwtVariogram->setGeometry(QRect(0, 10, 731, 391));

        retranslateUi(dlgSemivariogram);
        QObject::connect(buttonBox, SIGNAL(accepted()), dlgSemivariogram, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), dlgSemivariogram, SLOT(reject()));

        QMetaObject::connectSlotsByName(dlgSemivariogram);
    } // setupUi

    void setSamples(std::list<interp::InterpPoint> &samples) {
    	QwtSymbol *sym=new QwtSymbol(QwtSymbol::Ellipse, QBrush(Qt::white), QPen(Qt::black), QSize(5,5));
    	for(auto it = samples.begin(); it != samples.end(); ++it) {
    		QwtPlotMarker *m = new QwtPlotMarker("Test");
    		m->setSymbol(sym);
    		m->setValue(it->x, it->y);
    		m->attach(qwtVariogram);
    	}
    	qwtVariogram->replot();
    }

    void retranslateUi(QDialog *dlgSemivariogram)
    {
        dlgSemivariogram->setWindowTitle(QApplication::translate("dlgSemivariogram", "Dialog", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace interp {
namespace kriging {
namespace ui {
    class KrigePlot : public Ui_dlgSemivariogram {};
}
}
}

QT_END_NAMESPACE

#endif // VARIOGRAMYJ7917_H

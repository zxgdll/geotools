/********************************************************************************
** Form generated from reading UI file 'KrigePlot.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_KRIGEPLOT_H
#define UI_KRIGEPLOT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSlider>
#include <QtGui/QVBoxLayout>
#include "qwt_plot.h"

QT_BEGIN_NAMESPACE

class Ui_KrigePlot
{
public:
    QVBoxLayout *verticalLayout;
    QwtPlot *vPlot;
    QHBoxLayout *vFormLayout;
    QGridLayout *variogramFormGrid;
    QDoubleSpinBox *rangeDoubleSpinBox;
    QLabel *nuggetLabel;
    QLabel *vModelLabel;
    QLabel *sillLabel;
    QLabel *rangeLabel;
    QDoubleSpinBox *sillDoubleSpinBox;
    QDoubleSpinBox *nuggetDoubleSpinBox;
    QSlider *nuggetSlider;
    QSlider *sillSlider;
    QSlider *rangeSlider;
    QComboBox *vModelComboBox;
    QDialogButtonBox *okCancelGroup;

    void setupUi(QDialog *KrigePlot)
    {
        if (KrigePlot->objectName().isEmpty())
            KrigePlot->setObjectName(QString::fromUtf8("KrigePlot"));
        KrigePlot->resize(739, 587);
        KrigePlot->setModal(true);
        verticalLayout = new QVBoxLayout(KrigePlot);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        vPlot = new QwtPlot(KrigePlot);
        vPlot->setObjectName(QString::fromUtf8("vPlot"));
        vPlot->setAutoReplot(true);

        verticalLayout->addWidget(vPlot);

        vFormLayout = new QHBoxLayout();
        vFormLayout->setObjectName(QString::fromUtf8("vFormLayout"));
        variogramFormGrid = new QGridLayout();
        variogramFormGrid->setObjectName(QString::fromUtf8("variogramFormGrid"));
        rangeDoubleSpinBox = new QDoubleSpinBox(KrigePlot);
        rangeDoubleSpinBox->setObjectName(QString::fromUtf8("rangeDoubleSpinBox"));
        rangeDoubleSpinBox->setMaximum(1e+09);

        variogramFormGrid->addWidget(rangeDoubleSpinBox, 4, 3, 1, 1);

        nuggetLabel = new QLabel(KrigePlot);
        nuggetLabel->setObjectName(QString::fromUtf8("nuggetLabel"));

        variogramFormGrid->addWidget(nuggetLabel, 1, 0, 1, 1);

        vModelLabel = new QLabel(KrigePlot);
        vModelLabel->setObjectName(QString::fromUtf8("vModelLabel"));

        variogramFormGrid->addWidget(vModelLabel, 0, 0, 1, 1);

        sillLabel = new QLabel(KrigePlot);
        sillLabel->setObjectName(QString::fromUtf8("sillLabel"));

        variogramFormGrid->addWidget(sillLabel, 3, 0, 1, 1);

        rangeLabel = new QLabel(KrigePlot);
        rangeLabel->setObjectName(QString::fromUtf8("rangeLabel"));

        variogramFormGrid->addWidget(rangeLabel, 4, 0, 1, 1);

        sillDoubleSpinBox = new QDoubleSpinBox(KrigePlot);
        sillDoubleSpinBox->setObjectName(QString::fromUtf8("sillDoubleSpinBox"));
        sillDoubleSpinBox->setInputMethodHints(Qt::ImhFormattedNumbersOnly);
        sillDoubleSpinBox->setDecimals(3);
        sillDoubleSpinBox->setMaximum(1e+08);
        sillDoubleSpinBox->setValue(0);

        variogramFormGrid->addWidget(sillDoubleSpinBox, 3, 3, 1, 1);

        nuggetDoubleSpinBox = new QDoubleSpinBox(KrigePlot);
        nuggetDoubleSpinBox->setObjectName(QString::fromUtf8("nuggetDoubleSpinBox"));
        nuggetDoubleSpinBox->setDecimals(3);
        nuggetDoubleSpinBox->setMaximum(1e+08);

        variogramFormGrid->addWidget(nuggetDoubleSpinBox, 1, 3, 1, 1);

        nuggetSlider = new QSlider(KrigePlot);
        nuggetSlider->setObjectName(QString::fromUtf8("nuggetSlider"));
        nuggetSlider->setOrientation(Qt::Horizontal);

        variogramFormGrid->addWidget(nuggetSlider, 1, 2, 1, 1);

        sillSlider = new QSlider(KrigePlot);
        sillSlider->setObjectName(QString::fromUtf8("sillSlider"));
        sillSlider->setOrientation(Qt::Horizontal);

        variogramFormGrid->addWidget(sillSlider, 3, 2, 1, 1);

        rangeSlider = new QSlider(KrigePlot);
        rangeSlider->setObjectName(QString::fromUtf8("rangeSlider"));
        rangeSlider->setOrientation(Qt::Horizontal);

        variogramFormGrid->addWidget(rangeSlider, 4, 2, 1, 1);

        vModelComboBox = new QComboBox(KrigePlot);
        vModelComboBox->setObjectName(QString::fromUtf8("vModelComboBox"));

        variogramFormGrid->addWidget(vModelComboBox, 0, 2, 1, 2);


        vFormLayout->addLayout(variogramFormGrid);


        verticalLayout->addLayout(vFormLayout);

        okCancelGroup = new QDialogButtonBox(KrigePlot);
        okCancelGroup->setObjectName(QString::fromUtf8("okCancelGroup"));
        okCancelGroup->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(okCancelGroup);


        retranslateUi(KrigePlot);

        QMetaObject::connectSlotsByName(KrigePlot);
    } // setupUi

    void retranslateUi(QDialog *KrigePlot)
    {
        KrigePlot->setWindowTitle(QApplication::translate("KrigePlot", "Semivariogram Plot", 0, QApplication::UnicodeUTF8));
        nuggetLabel->setText(QApplication::translate("KrigePlot", "Nugget", 0, QApplication::UnicodeUTF8));
        vModelLabel->setText(QApplication::translate("KrigePlot", "Variogram Model", 0, QApplication::UnicodeUTF8));
        sillLabel->setText(QApplication::translate("KrigePlot", "Sill", 0, QApplication::UnicodeUTF8));
        rangeLabel->setText(QApplication::translate("KrigePlot", "Range", 0, QApplication::UnicodeUTF8));
        vModelComboBox->clear();
        vModelComboBox->insertItems(0, QStringList()
         << QApplication::translate("KrigePlot", "Spherical", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("KrigePlot", "Gaussian", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("KrigePlot", "Linear", 0, QApplication::UnicodeUTF8)
        );
    } // retranslateUi

};

namespace Ui {
    class KrigePlot: public Ui_KrigePlot {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_KRIGEPLOT_H

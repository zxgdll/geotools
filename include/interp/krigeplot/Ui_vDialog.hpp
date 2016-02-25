/********************************************************************************
** Form generated from reading UI file 'variogramf11779.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef VARIOGRAMF11779_H
#define VARIOGRAMF11779_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFormLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QVBoxLayout>
#include "qwt_plot.h"

QT_BEGIN_NAMESPACE

class Ui_vDialog
{
public:
    QVBoxLayout *verticalLayout;
    QwtPlot *vPlot;
    QHBoxLayout *vFormLayout;
    QFormLayout *vLeftFormLayout;
    QLabel *vModelLabel;
    QComboBox *vModelComboBox;
    QLabel *nuggetLabel;
    QLabel *sillLabel;
    QDoubleSpinBox *nuggetDoubleSpinBox;
    QDoubleSpinBox *sillDoubleSpinBox;
    QFormLayout *vRightFormLayout;
    QDialogButtonBox *okCancelGroup;

    void setupUi(QDialog *vDialog)
    {
        if (vDialog->objectName().isEmpty())
            vDialog->setObjectName(QString::fromUtf8("vDialog"));
        vDialog->resize(739, 587);
        vDialog->setModal(true);
        verticalLayout = new QVBoxLayout(vDialog);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        vPlot = new QwtPlot(vDialog);
        vPlot->setObjectName(QString::fromUtf8("vPlot"));
        vPlot->setAutoReplot(true);

        verticalLayout->addWidget(vPlot);

        vFormLayout = new QHBoxLayout();
        vFormLayout->setObjectName(QString::fromUtf8("vFormLayout"));
        vLeftFormLayout = new QFormLayout();
        vLeftFormLayout->setObjectName(QString::fromUtf8("vLeftFormLayout"));
        vLeftFormLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
        vModelLabel = new QLabel(vDialog);
        vModelLabel->setObjectName(QString::fromUtf8("vModelLabel"));

        vLeftFormLayout->setWidget(0, QFormLayout::LabelRole, vModelLabel);

        vModelComboBox = new QComboBox(vDialog);
        vModelComboBox->setObjectName(QString::fromUtf8("vModelComboBox"));

        vLeftFormLayout->setWidget(0, QFormLayout::FieldRole, vModelComboBox);

        nuggetLabel = new QLabel(vDialog);
        nuggetLabel->setObjectName(QString::fromUtf8("nuggetLabel"));

        vLeftFormLayout->setWidget(1, QFormLayout::LabelRole, nuggetLabel);

        sillLabel = new QLabel(vDialog);
        sillLabel->setObjectName(QString::fromUtf8("sillLabel"));

        vLeftFormLayout->setWidget(2, QFormLayout::LabelRole, sillLabel);

        nuggetDoubleSpinBox = new QDoubleSpinBox(vDialog);
        nuggetDoubleSpinBox->setObjectName(QString::fromUtf8("nuggetDoubleSpinBox"));
        nuggetDoubleSpinBox->setDecimals(3);
        nuggetDoubleSpinBox->setMaximum(1e+08);

        vLeftFormLayout->setWidget(1, QFormLayout::FieldRole, nuggetDoubleSpinBox);

        sillDoubleSpinBox = new QDoubleSpinBox(vDialog);
        sillDoubleSpinBox->setObjectName(QString::fromUtf8("sillDoubleSpinBox"));
        sillDoubleSpinBox->setInputMethodHints(Qt::ImhFormattedNumbersOnly);
        sillDoubleSpinBox->setDecimals(3);
        sillDoubleSpinBox->setMaximum(1e+08);
        sillDoubleSpinBox->setValue(1);

        vLeftFormLayout->setWidget(2, QFormLayout::FieldRole, sillDoubleSpinBox);


        vFormLayout->addLayout(vLeftFormLayout);

        vRightFormLayout = new QFormLayout();
        vRightFormLayout->setObjectName(QString::fromUtf8("vRightFormLayout"));

        vFormLayout->addLayout(vRightFormLayout);


        verticalLayout->addLayout(vFormLayout);

        okCancelGroup = new QDialogButtonBox(vDialog);
        okCancelGroup->setObjectName(QString::fromUtf8("okCancelGroup"));
        okCancelGroup->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(okCancelGroup);


        retranslateUi(vDialog);

        QMetaObject::connectSlotsByName(vDialog);
    } // setupUi

    void retranslateUi(QDialog *vDialog)
    {
        vDialog->setWindowTitle(QApplication::translate("vDialog", "Semivariogram Plot", 0, QApplication::UnicodeUTF8));
        vModelLabel->setText(QApplication::translate("vDialog", "Variogram Model", 0, QApplication::UnicodeUTF8));
        vModelComboBox->clear();
        vModelComboBox->insertItems(0, QStringList()
         << QApplication::translate("vDialog", "Spherical", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("vDialog", "Gaussian", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("vDialog", "Linear", 0, QApplication::UnicodeUTF8)
        );
        nuggetLabel->setText(QApplication::translate("vDialog", "Nugget", 0, QApplication::UnicodeUTF8));
        sillLabel->setText(QApplication::translate("vDialog", "Sill", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class vDialog: public Ui_vDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // VARIOGRAMF11779_H

/********************************************************************************
** Form generated from reading UI file 'variogramyJ7917.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef _KRIGEPLOT_HPP_
#define _KRIGEPLOT_HPP_

#include "qwt_plot.h"
#include "qwt_plot_marker.h"
#include "qwt_symbol.h"

#include <limits>
#include <list>

#include "interp/krigeplot/Ui_vDialog.hpp"
#include "interp/SimpleKrigingInterpolator.hpp"
#include "interp/InterpPoint.hpp"

class _KrigePlot : public Ui_vDialog {

public:
	void setVariogram(std::list<interp::kriging::VariogramPoint> &variogram) {
		double maxx, maxy;
		maxx = maxy = std::numeric_limits<double>::lowest();
		for(auto it = variogram.begin(); it != variogram.end(); ++it) {
			QwtSymbol *sym=new QwtSymbol(QwtSymbol::Ellipse, QBrush(Qt::white), QPen(Qt::black), QSize(5,5));
			QwtPlotMarker *m = new QwtPlotMarker("Test");
			m->setSymbol(sym);
			m->setValue(it->distance(), it->difference());
			m->attach(vPlot);
			if(it->distance() > maxx) maxx = it->distance();
			if(it->difference() > maxy) maxy = it->difference();
			//std::cerr << it->distance() << ", " << it->difference() << std::endl;
		}
		vPlot->setAxisScale(QwtPlot::xBottom, 0, maxx);
		vPlot->setAxisScale(QwtPlot::yLeft, 0, maxy);
		vPlot->replot();
	}

	void vModelIndexChanged(int index) {
		std::cerr << "combo " << index << std::endl;
	}

	void setupUi(QDialog *vDialog) {
		Ui_vDialog::setupUi(vDialog);
		//QObject::connect(vModelComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(vModelIndexChanged(int)));
	}

};

namespace interp {
	namespace kriging {
		namespace ui {

			class KrigePlot : public _KrigePlot {
			};

		}
	}
}

#endif

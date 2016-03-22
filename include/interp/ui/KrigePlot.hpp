#ifndef _KRIGEPLOT_HPP_
#define _KRIGEPLOT_HPP_

#include <QtCore/QObject>
#include "qwt_plot.h"
#include "qwt_plot_marker.h"
#include "qwt_plot_curve.h"
#include "qwt_symbol.h"

#include <iostream>
#include <limits>
#include <list>

#include "interp/ui/ui_KrigePlot.h"
#include "interp/SimpleKrigingInterpolator.hpp"
#include "interp/InterpPoint.hpp"

namespace interp {

	namespace ui {

		class KrigePlot :  public QObject, public Ui::KrigePlot {

		Q_OBJECT

		private:

			// O'Sullivan, Unwin p. 298.
			static double spherical(double x, double nug, double sill, double range) {
				if(x > range) {
					return nug + sill;
				} else {
					return nug + sill * ((3 * x) / (2 * range) - 0.5 * pow((x / range), 3.0));
				}
			}

			static double linear(double x, double nug, double sill, double range) {
				return (sill - nug) / range * x + nug;
			}

			static double gaussian(double x, double nug, double sill, double range) {
				return linear(x, nug, sill, range);
			}

			double (* m_models[3])(double, double, double, double);
			double (* m_model)(double, double, double, double);

			interp::kriging::KrigeArgs *m_args;

			double m_maxx;
			double m_maxy;

			QwtPlotCurve *m_curve;
			QDialog *m_vDialog;

			std::list<QwtPlotItem*> m_items;
			std::list<interp::kriging::VariogramPoint> m_variogram;

			KrigePlot() : Ui::KrigePlot() {
				m_models[0] = &(this->spherical);
				m_models[1] = &(this->gaussian);
				m_models[2] = &(this->linear);
				m_model = m_models[0];
				m_curve = nullptr;
				m_maxx = 0.0;
				m_maxy = 0.0;
				m_args = nullptr;
				m_vDialog = nullptr;
			}

			void drawVariogram() {
				std::cerr << m_variogram.size() << " points" << std::endl;
				if(m_items.size() == 0) {
					m_maxx = std::numeric_limits<double>::lowest();
					m_maxy = std::numeric_limits<double>::lowest();
					for(auto it = m_variogram.begin(); it != m_variogram.end(); ++it) {
						QwtSymbol *sym=new QwtSymbol(QwtSymbol::Ellipse, QBrush(Qt::white), QPen(Qt::black), QSize(5,5));
						QwtPlotMarker *m = new QwtPlotMarker("Test");
						m->setSymbol(sym);
						m->setValue(it->distance(), it->difference());
						m->setZ(0.0);
						m->attach(vPlot);
						m_items.push_back(m);
						//std::cerr << it->distance() << std::endl;
						if(it->distance() > m_maxx)
							m_maxx = it->distance();
						if(it->difference() > m_maxy)
							m_maxy = it->difference();
					}

					m_args->sill = m_maxy / 2.0;
					m_args->range = m_maxx / 2.0;
					m_args->model = m_model;

					std::cerr << m_args->sill << ", " << m_args->range << std::endl;

					vPlot->setAxisScale(QwtPlot::xBottom, 0, m_maxx);
					vPlot->setAxisScale(QwtPlot::yLeft, 0, m_maxy);

					nuggetDoubleSpinBox->setValue(m_args->nugget);
					nuggetSlider->setRange(0, m_maxy * 100);
					nuggetSlider->setValue(m_args->nugget * 100);

					sillDoubleSpinBox->setValue(m_args->sill);
					sillSlider->setRange(0, m_maxy * 100);
					sillSlider->setValue(m_args->sill * 100);

					rangeDoubleSpinBox->setValue(m_args->range);
					rangeSlider->setRange(0, m_maxx * 100);
					rangeSlider->setValue(m_args->range * 100);
				}

				if(m_curve == nullptr) {
					m_curve = new QwtPlotCurve("Model");
					m_curve->setPen(QPen(Qt::red));
					m_curve->setZ(1.0);
					m_curve->attach(vPlot);
				}

				std::vector<double> xSamples(100);
				std::vector<double> ySamples(100);
				for(double x = 0; x < m_maxx; x += m_maxx / 100.0) {
					xSamples.push_back(x);
					ySamples.push_back((double) m_model(x, m_args->nugget, m_args->sill, m_args->range));
				}
				m_curve->setSamples(xSamples.data(), ySamples.data(), (int) xSamples.size());

				vPlot->replot();
			}

			void valueUpdate() {
				m_args->nugget = nuggetDoubleSpinBox->value();
				m_args->sill = sillDoubleSpinBox->value();
				m_args->range = rangeDoubleSpinBox->value();
				drawVariogram();
			}


		public slots:

			void vModelIndexChanged(int index) {
				m_model = m_models[index];
				valueUpdate();
			}

			void sillChanged() {
				sillSlider->setValue(sillDoubleSpinBox->value() * 100);
				valueUpdate();
			}

			void sillSliderChanged(int value) {
				sillDoubleSpinBox->setValue(value / 100.0);
				valueUpdate();
			}

			void rangeSliderChanged(int value) {
				rangeDoubleSpinBox->setValue(value / 100.0);
				valueUpdate();
			}

			void rangeChanged() {
				rangeSlider->setValue(rangeDoubleSpinBox->value() * 100);
				valueUpdate();
			}

			void nuggetSliderChanged(int value) {
				nuggetDoubleSpinBox->setValue(value / 100.0);
				valueUpdate();
			}

			void nuggetChanged() {
				nuggetSlider->setValue(nuggetDoubleSpinBox->value() * 100);
				valueUpdate();
			}

			void cancel() {
				m_args->status = 1;
				m_vDialog->done(1);
			}

			void ok() {
				m_args->status = 0;
				m_vDialog->done(0);
			}

		public:

			KrigePlot(interp::kriging::KrigeArgs *args) : KrigePlot() {
				m_args = args;
			}

			void setVariogram(std::list<interp::kriging::VariogramPoint> &variogram) {
				m_variogram.clear();
				m_variogram.assign(variogram.begin(), variogram.end());
				drawVariogram();
			}

			void setupUi(QDialog *vDialog) {
				m_vDialog = vDialog;
				Ui::KrigePlot::setupUi(vDialog);
				QObject::connect(vModelComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(vModelIndexChanged(int)));

				QObject::connect(sillSlider, SIGNAL(valueChanged(int)), this, SLOT(sillSliderChanged(int)));
				QObject::connect(sillDoubleSpinBox, SIGNAL(editingFinished()), this, SLOT(sillChanged()));

				QObject::connect(rangeSlider, SIGNAL(valueChanged(int)), this, SLOT(rangeSliderChanged(int)));
				QObject::connect(rangeDoubleSpinBox, SIGNAL(editingFinished()), this, SLOT(rangeChanged()));

				QObject::connect(nuggetSlider, SIGNAL(valueChanged(int)), this, SLOT(nuggetSliderChanged(int)));
				QObject::connect(nuggetDoubleSpinBox, SIGNAL(editingFinished()), this, SLOT(nuggetChanged()));

				QObject::connect(okCancelGroup, SIGNAL(accepted()), this, SLOT(ok()));
				QObject::connect(okCancelGroup, SIGNAL(rejected()), this, SLOT(cancel()));
			}
		};


	} // ui
} // interp

#endif

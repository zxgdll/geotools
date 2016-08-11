#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <QtWidgets/QWidget>

#include "ui_lasgrid.h"

namespace geotools {

	namespace ui {

		class LasgridForm : public QObject, public Ui::LasgridForm {
			Q_OBJECT
		private:
			QWidget *m_form;
		public:
			LasgridForm(QWidget *p = 0);
			void setupUi(QWidget *Form);
		public slots:
			void selectFileClicked();
		};

	}

}

#endif


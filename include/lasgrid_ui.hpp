#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <QtGui/QWidget>

#include "ui_lasgrid.h"

namespace geotools {

	namespace ui {

		class LasgridForm : public Ui::LasgridForm {
		public:
			void setupUi(QWidget *Form);
		};

	}

}

#endif


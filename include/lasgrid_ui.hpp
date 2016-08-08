#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <QtGui/QWidget>

#include "ui_lasgrid.h"

namespace geotools {

	namespace ui {

		class LasgridForm : public Ui::Form {
		public:
			void setupUi(QWidget *Form) {
	                        Ui::Form::setupUi(Form);
        	                int i = 0;
                	        std::map<std::string, int> types = geotools::las::config::types;
                        	for(auto it = types.begin(); it != types.end(); ++it)
                                	cboType->insertItem(i++, QString::fromStdString(it->first), QVariant(it->second));
			}
		};

	}

}

#endif


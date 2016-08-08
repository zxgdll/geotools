#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <QtGui/QWidget>

namespace geotools {

	namespace ui {

		class LasgridForm : public Ui::Form {
		public:
			void setupUi(QWidget *Form) {
				using namespace geotools::las::config;

	                        Ui::Form::setupUi(Form);

        	                int i = 0;
				int defaultIdx = -1;
                        	for(auto it = types.begin(); it != types.end(); ++it) {
                                	cboType->insertItem(i, QString::fromStdString(it->first), QVariant(it->second));
					if(it->second == defaultType)
						defaultIdx = i;
					++i;
				}
				cboType->setCurrentIndex(defaultIdx);

				i = 0;
				defaultIdx = -1;
				std::map<std::string, int> atts = geotools::las::config::attributes;
				for(auto it = atts.begin(); it != atts.end(); ++it) {
					cboAttribute->insertItem(i, QString::fromStdString(it->first), QVariant(it->second));
					if(it->second == defaultAttribute)
						defaultIdx = i;
					++i;
				}
				cboAttribute->setCurrentIndex(defaultIdx);
			}
		};

	}

}

#endif


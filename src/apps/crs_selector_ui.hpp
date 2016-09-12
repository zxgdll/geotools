#ifndef __CRS_SELECTOR_UI_HPP__
#define __CRS_SELECTOR_UI_HPP__

#include <map>
#include <string>

#include <QtWidgets/QWidget>
#include <QDir>

#include "ui_crs_selector.h"

namespace geotools {

        namespace ui {

                class CRSSelector : public QDialog, public Ui::CRSSelector {
                        Q_OBJECT
                private:
			std::map<int, std::string> m_hcrs;
			std::map<int, std::string> m_vcrs;
                        void initUi();
			void loadCrs(std::map<int, std::string> &target, const std::string &filename);
                public:
                       	CRSSelector(QWidget *parent = Q_NULLPTR, Qt::WindowFlags f = Qt::WindowFlags());

                public slots:
			void test();
                };

        }

}

#endif



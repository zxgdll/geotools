#ifndef __LASGRID_UI_HPP__
#define __LASGRID_UI_HPP__

#include <QtWidgets/QWidget>
#include <QDir>

#include "ui_lasgrid.h"

namespace geotools {

	namespace ui {

		class LasgridForm : public QWidget, public Ui::LasgridForm {
			Q_OBJECT
		private:
			QWidget *m_form;
			std::string m_destFile;
			std::vector<std::string> m_lasFiles;
			std::set<int> m_classes;
			int m_crs;
			int m_attribute;
			int m_type;
			unsigned char m_angleLimit;
			bool m_fill;
			bool m_snap;
			double m_resolution;
			double m_radius;
			QDir *m_last;
			void updateFileList();
			void updateFileButtons();
		public:
			LasgridForm(QWidget *p = Q_NULLPTR);
			void setupUi(QWidget *Form);
			~LasgridForm();
		public slots:
			void fileListSelectionChanged();
			void selectFilesClicked();
			void removeFilesClicked();
			void clearFilesClicked();
			void cancelClicked();
			void runClicked();
			void destFileClicked();
			void crsConfigClicked();
		};

	}

}

#endif


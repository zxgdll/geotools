#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>

#include "csv.h"

#include "geotools.h"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;

CRSSelector::CRSSelector(QWidget *parent, Qt::WindowFlags f) : 
	QDialog(parent, f) {
	initUi();
}

void CRSSelector::loadCrs(std::map<int, std::string> &target, const std::string &filename) {
        io::CSVReader<2, io::trim_chars<' ','\t'>, io::double_quote_escape<',','"'>, io::throw_on_overflow, io::single_line_comment<'#'> > in(filename.c_str());
        in.read_header(io::ignore_extra_column, "COORD_REF_SYS_CODE","COORD_REF_SYS_NAME");
        int code;
	std::string name;
        while(in.read_row(code, name))
		target[code] = name;
}

void CRSSelector::initUi() {
	Ui:CRSSelector::setupUi((QDialog *) this);

	std::map<int, std::string> vcrs;
	loadCrs(vcrs, "/usr/share/gdal/1.11/vertcs.csv");
	loadCrs(vcrs, "/usr/share/gdal/1.11/vertcs.override.csv");

	for(const std::pair<int, std::string> &item : vcrs)
		cboVerticalDatum->addItem(QString(item.second.c_str()), QVariant(item.first));

}

void CRSSelector::test() {}



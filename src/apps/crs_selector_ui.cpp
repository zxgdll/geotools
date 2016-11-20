#include <iterator>

#include <QtWidgets/QWidget>
#include <QtWidgets/QFileDialog>
#include <QCompleter>

#include "cpl_conv.h"

#include "csv.h"

#include "geotools.h"
#include "crs_selector_ui.hpp"

using namespace geotools::ui;

CRSSelector::CRSSelector(QWidget *parent, Qt::WindowFlags f) :
QDialog(parent, f),
m_vsrid(0), m_hsrid(0) {
    initUi();
}

void CRSSelector::loadCrs(std::map<int, std::string> &target, const std::string &filename) {
    const char *path = CPLFindFile("gdal", filename.c_str());
    if (path == NULL)
        g_argerr("The path to " << filename << " could not be determined. Please set GDAL_DATA.");
    io::CSVReader<2, io::trim_chars<' ', '\t'>, io::double_quote_escape<',', '"'>, io::throw_on_overflow, io::single_line_comment<'#'> > in(path);
    in.read_header(io::ignore_extra_column, "COORD_REF_SYS_CODE", "COORD_REF_SYS_NAME");
    int code;
    std::string name;
    while (in.read_row(code, name))
        target[code] = name;
}

void CRSSelector::initUi() {
    Ui::CRSSelector::setupUi((QDialog *) this);

    loadCrs(m_vcrs, "vertcs.csv");
    loadCrs(m_vcrs, "vertcs.override.csv");
    for (const std::pair<int, std::string> &item : m_vcrs)
        m_vlst << QString((std::to_string(item.first) + " - " + item.second).c_str());

    m_vcomp = new QCompleter(m_vlst, this);
    m_vcomp->setCaseSensitivity(Qt::CaseInsensitive);
    m_vcomp->setFilterMode(Qt::MatchContains);
    txtVerticalDatum->setCompleter(m_vcomp);

    loadCrs(m_hcrs, "gcs.csv");
    loadCrs(m_hcrs, "pcs.csv");
    for (const std::pair<int, std::string> &item : m_hcrs)
        m_hlst << QString((std::to_string(item.first) + " - " + item.second).c_str());

    m_hcomp = new QCompleter(m_hlst, this);
    m_hcomp->setCaseSensitivity(Qt::CaseInsensitive);
    m_hcomp->setFilterMode(Qt::MatchContains);
    txtHorizontalCRS->setCompleter(m_hcomp);

    connect(btnSelect, SIGNAL(clicked()), this, SLOT(selectClicked()));
    connect(btnCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));
    connect(txtHorizontalCRS, SIGNAL(textChanged(const QString &)), this, SLOT(textUpdate()));
    connect(txtVerticalDatum, SIGNAL(textChanged(const QString &)), this, SLOT(textUpdate()));

    updateFields();
}

void CRSSelector::updateFields() {
    if (m_hcomp && m_hsrid > 0) {
        QString h;
        h.setNum(m_hsrid);
        txtHorizontalCRS->setText(h);
        m_hcomp->complete();
    }
    if (m_vcomp && m_vsrid > 0) {
        QString v;
        v.setNum(m_vsrid);
        txtVerticalDatum->setText(v);
        m_vcomp->complete();
    }
}

void CRSSelector::textUpdate() {
    btnSelect->setEnabled(m_hcomp->currentIndex().isValid());
}

void CRSSelector::cancelClicked() {
    reject();
}

void CRSSelector::selectClicked() {
    m_hsrid = 0;
    m_vsrid = 0;
    int hidx = m_hlst.indexOf(m_hcomp->currentCompletion());
    if (hidx >= 0)
        m_hsrid = std::next(m_hcrs.begin(), hidx)->first;
    int vidx = m_vlst.indexOf(m_vcomp->currentCompletion());
    if (vidx >= 0)
        m_vsrid = std::next(m_vcrs.begin(), vidx)->first;
    g_debug("vsrid: " << m_vsrid << "; hsrid: " << m_hsrid);
    accept();
}

void CRSSelector::setVerticalSRID(int srid) {
    m_vsrid = srid;
    updateFields();
}

void CRSSelector::setHorizontalSRID(int srid) {
    m_hsrid = srid;
    updateFields();
}

int CRSSelector::getVerticalSRID() {
    return m_vsrid;
}

int CRSSelector::getHorizontalSRID() {
    return m_hsrid;
}




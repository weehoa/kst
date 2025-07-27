/***************************************************************************
				netcdf.cpp  -  netCDF file data source reader
							 -------------------
	begin                : 17/06/2004
	copyright            : (C) 2004 Nicolas Brisset <nicodev@users.sourceforge.net>
	email                : kst@kde.org
	modified             : 03/14/05 by K. Scott
 ***************************************************************************/

 /***************************************************************************
  *                                                                         *
  *   This program is free software; you can redistribute it and/or modify  *
  *   it under the terms of the GNU General Public License as published by  *
  *   the Free Software Foundation; either version 2 of the License, or     *
  *   (at your option) any later version.                                   *
  *                                                                         *
  ***************************************************************************/

#include "sharedptr.h"
#include "netcdfsource.h"
#include "debug.h"

#include <QFile>
#include <QFileInfo>
#include <ctype.h>
#include <stdlib.h>
#include <netcdf>
#include <memory>
#include <map>
#include <vector>
#include <string>

using namespace netCDF;
using namespace Kst;

static const QString netCdfTypeString = "netCDF Files";
static constexpr const char* PATH_SEP = "/";
static constexpr const char* INDEX_FIELD = "INDEX";
static constexpr const char* INDEX_SUFFIX = "/INDEX";


//
// Scalar interface
//

class DataInterfaceNetCdfScalar : public DataSource::DataInterface<DataScalar>
{
public:
	DataInterfaceNetCdfScalar(NetcdfSource& s) : netcdf(s) {}

	// read one element
	int read(const QString&, DataScalar::ReadInfo&);

	// named elements
	QStringList list() const { return netcdf._scalarList; }
	bool isListComplete() const { return true; }
	bool isValid(const QString&) const;

	// T specific
	const DataScalar::DataInfo dataInfo(const QString&, int frame = 0) const { Q_UNUSED(frame) return DataScalar::DataInfo(); }
	void setDataInfo(const QString&, const DataScalar::DataInfo&) {}

	// meta data
	QMap<QString, double> metaScalars(const QString&) { return QMap<QString, double>(); }
	QMap<QString, QString> metaStrings(const QString&) { return QMap<QString, QString>(); }


private:
	NetcdfSource& netcdf;
};

int DataInterfaceNetCdfScalar::read(const QString& scalar, DataScalar::ReadInfo& p)
{
	return netcdf.readScalar(p.value, scalar);
}

bool DataInterfaceNetCdfScalar::isValid(const QString& scalar) const
{
	return  netcdf._scalarList.contains(scalar);
}

//
// String interface
//

class DataInterfaceNetCdfString : public DataSource::DataInterface<DataString>
{
public:
	DataInterfaceNetCdfString(NetcdfSource& s) : netcdf(s) {}

	// read one element
	int read(const QString&, DataString::ReadInfo&);

	// named elements
	QStringList list() const { return netcdf._strings.keys(); }
	bool isListComplete() const { return true; }
	bool isValid(const QString&) const;

	// T specific
	const DataString::DataInfo dataInfo(const QString&, int frame = 0) const { Q_UNUSED(frame) return DataString::DataInfo(); }
	void setDataInfo(const QString&, const DataString::DataInfo&) {}

	// meta data
	QMap<QString, double> metaScalars(const QString&) { return QMap<QString, double>(); }
	QMap<QString, QString> metaStrings(const QString&) { return QMap<QString, QString>(); }


private:
	NetcdfSource& netcdf;
};

int DataInterfaceNetCdfString::read(const QString& string, DataString::ReadInfo& p)
{
	if (isValid(string) && p.value) {
		*p.value = netcdf._strings[string];
		return 1;
	}
	return 0;
}

bool DataInterfaceNetCdfString::isValid(const QString& string) const
{
	return netcdf._strings.contains(string);
}

//
// Vector interface
//

class DataInterfaceNetCdfVector : public DataSource::DataInterface<DataVector>
{
public:
	DataInterfaceNetCdfVector(NetcdfSource& s) : netcdf(s) {}

	int read(const QString&, DataVector::ReadInfo&);
	QStringList list() const { return netcdf._fieldList; }
	bool isListComplete() const { return true; }
	bool isValid(const QString&) const;
	const DataVector::DataInfo dataInfo(const QString&, int frame = 0) const;
	void setDataInfo(const QString&, const DataVector::DataInfo&) {}
	QMap<QString, double> metaScalars(const QString&);
	QMap<QString, QString> metaStrings(const QString&);

private:
	NetcdfSource& netcdf;
};

const DataVector::DataInfo DataInterfaceNetCdfVector::dataInfo(const QString& field, int frame) const
{
	Q_UNUSED(frame)
	if (!netcdf._fieldList.contains(field))
		return DataVector::DataInfo();

	return DataVector::DataInfo(netcdf.frameCount(field), netcdf.samplesPerFrame(field));
}

int DataInterfaceNetCdfVector::read(const QString& field, DataVector::ReadInfo& p)
{
	return netcdf.readField(p.data, field, p.startingFrame, p.numberOfFrames);
}

bool DataInterfaceNetCdfVector::isValid(const QString& field) const
{
	return netcdf._fieldList.contains(field);
}

QMap<QString, double> DataInterfaceNetCdfVector::metaScalars(const QString& field)
{
	QMap<QString, double> fieldScalars;
	if (!netcdf._ncfile || NetcdfSource::isIndexField(field)) {
		return fieldScalars;
	}
	try {
		NcVar var = netcdf.getVarByFullPath(field);
		if (var.isNull()) {
			qDebug() << "Queried field " << field << " which can't be read" << endl;
			return fieldScalars;
		}
		auto atts = var.getAtts();
		fieldScalars["NbAttributes"] = static_cast<double>(atts.size());
		for (const auto& attPair : atts) {
			const NcAtt& att = attPair.second;
			NcType attType = att.getType();
			NcType::ncType t = attType.getTypeClass();
			if (t == NC_BYTE || t == NC_SHORT ||
				t == NC_INT || t == NC_INT64 ||
				t == NC_FLOAT || t == NC_DOUBLE) {
				std::vector<double> values;
				values.resize(att.getAttLength());
				att.getValues(values.data());
				if (!values.empty()) {
					fieldScalars[QString::fromStdString(att.getName())] = values[0];
					for (size_t j = 1; j < values.size(); ++j) {
						fieldScalars[QString::fromStdString(att.getName()) + QString("-%1").arg(j + 1)] = values[j];
					}
				}
			}
		}
	}
	catch (...) {}
	return fieldScalars;
}

QMap<QString, QString> DataInterfaceNetCdfVector::metaStrings(const QString& field)
{
	QMap<QString, QString> fieldStrings;
	if (!netcdf._ncfile || NetcdfSource::isIndexField(field)) {
		return fieldStrings;
	}
	try {
		NcVar var = netcdf.getVarByFullPath(field);
		if (var.isNull()) {
			qDebug() << "Queried field " << field << " which can't be read" << endl;
			return fieldStrings;
		}
		auto atts = var.getAtts();
		for (const auto& attPair : atts) {
			const NcAtt& att = attPair.second;
			NcType attType = att.getType();
			NcType::ncType t = attType.getTypeClass();
			if (t == NC_CHAR || t == NC_STRING) {
				std::string strVal;
				try {
					att.getValues(strVal);
					fieldStrings[QString::fromStdString(att.getName())] = QString::fromStdString(strVal);
				}
				catch (...) {}
			}
		}
	}
	catch (...) {}
	return fieldStrings;
}

//
// Matrix interface
//

class DataInterfaceNetCdfMatrix : public DataSource::DataInterface<DataMatrix>
{
public:
	DataInterfaceNetCdfMatrix(NetcdfSource& s) : netcdf(s) {}

	int read(const QString&, DataMatrix::ReadInfo&);
	QStringList list() const { return netcdf._matrixList; }
	bool isListComplete() const { return true; }
	bool isValid(const QString&) const;
	const DataMatrix::DataInfo dataInfo(const QString&, int frame) const;
	void setDataInfo(const QString&, const DataMatrix::DataInfo&) {}
	QMap<QString, double> metaScalars(const QString&) { return QMap<QString, double>(); }
	QMap<QString, QString> metaStrings(const QString&) { return QMap<QString, QString>(); }

private:
	NetcdfSource& netcdf;
};

const DataMatrix::DataInfo DataInterfaceNetCdfMatrix::dataInfo(const QString& matrix, int frame) const
{
	Q_UNUSED(frame)
	if (netcdf._ncfile || !netcdf._matrixList.contains(matrix)) {
		return DataMatrix::DataInfo();
	}
	try {
		NcVar var = netcdf.getVarByFullPath(matrix);
		if (var.isNull() || var.getDimCount() != 2) {
			return DataMatrix::DataInfo();
		}
		DataMatrix::DataInfo info;
		info.xSize = var.getDim(0).getSize();
		info.ySize = var.getDim(1).getSize();
		return info;
	}
	catch (...) {
		return DataMatrix::DataInfo();
	}
}

int DataInterfaceNetCdfMatrix::read(const QString& field, DataMatrix::ReadInfo& p)
{
	int count = netcdf.readMatrix(p.data->z, field);
	p.data->xMin = 0;
	p.data->yMin = 0;
	p.data->xStepSize = 1;
	p.data->yStepSize = 1;
	return count;
}

bool DataInterfaceNetCdfMatrix::isValid(const QString& field) const {
	return netcdf._matrixList.contains(field);
}

//
// NetcdfSource
//

NetcdfSource::NetcdfSource(Kst::ObjectStore* store, QSettings* cfg, const QString& filename, const QString& type, const QDomElement& element) :
	Kst::DataSource(store, cfg, filename, type),
	_ncfile(nullptr),
	is(new DataInterfaceNetCdfScalar(*this)),
	it(new DataInterfaceNetCdfString(*this)),
	iv(new DataInterfaceNetCdfVector(*this)),
	im(new DataInterfaceNetCdfMatrix(*this))
{
	setInterface(is);
	setInterface(it);
	setInterface(iv);
	setInterface(im);

	setUpdateType(None);

	if (!type.isEmpty() && type != "netCDF") {
		return;
	}

	_valid = false;

	_filename = filename;
	_strings = fileMetas();
	_valid = initFile();
}

NetcdfSource::~NetcdfSource() {
	_ncfile.reset();
}

void NetcdfSource::reset() {
	_ncfile.reset();
	_groupMaxFrameCount.clear();
	_valid = initFile();
}

// Returns the NcVar for a full path (with '/' separator), or a null NcVar if not found.
// Now uses the cached variable map.
NcVar NetcdfSource::getVarByFullPath(const QString& fullPath) const {
	auto it = _varMap.find(fullPath);
	if (it != _varMap.end()) {
		return it.value();
	}
	return NcVar(); // null NcVar
}

// In traverseGroups, compute and store per-group max frame count:
void NetcdfSource::traverseGroups(const NcGroup& group, const QString& path) {
	int groupMax = 0;
	auto vars = group.getVars();

	// Add an INDEX field for the current group
	QString indexField = path.isEmpty() ? INDEX_FIELD : path + INDEX_SUFFIX;
	_fieldList += indexField;
	_varMap[indexField] = NcVar();

	for (const auto& varPair : vars) {
		const NcVar& var = varPair.second;
		int numDims = var.getDimCount();
		QString fullPath = path.isEmpty() ? QString::fromStdString(var.getName())
			: path + PATH_SEP + QString::fromStdString(var.getName());

		if (numDims == 0) {
			_scalarList += fullPath;
		}
		else if (numDims <= 2) {
			size_t fc = var.getDim(0).getSize();
			if (numDims == 1) {
				_fieldList += fullPath;
			}
			else if (numDims == 2) {
				_matrixList += fullPath;
			}
			_frameCounts[fullPath] = static_cast<int>(fc);
			groupMax = qMax(groupMax, static_cast<int>(fc));
		}
		_varMap[fullPath] = var;
	}

	// Store per-group max frame count
	_groupMaxFrameCount[path] = groupMax;

	// Recursively traverse child groups
	auto groups = group.getGroups();
	for (const auto& groupPair : groups) {
		const NcGroup& childGroup = groupPair.second;
		QString childPath = path.isEmpty() ? QString::fromStdString(groupPair.first)
			: path + PATH_SEP + QString::fromStdString(groupPair.first);
		traverseGroups(childGroup, childPath);
	}
}

// Recursively traverse all groups and cache variables in _varMap.
void NetcdfSource::traverseAndCacheVars(const NcGroup& group, const QString& path) {
	auto vars = group.getVars();
	for (const auto& varPair : vars) {
		const NcVar& var = varPair.second;
		QString fullPath = path.isEmpty() ? QString::fromStdString(var.getName())
			: path + PATH_SEP + QString::fromStdString(var.getName());
		_varMap[fullPath] = var;
	}
	auto groups = group.getGroups();
	for (const auto& groupPair : groups) {
		const NcGroup& childGroup = groupPair.second;
		QString childPath = path.isEmpty() ? QString::fromStdString(groupPair.first)
			: path + PATH_SEP + QString::fromStdString(groupPair.first);
		traverseAndCacheVars(childGroup, childPath);
	}
}

void NetcdfSource::traverseAndUpdateGroups(const NcGroup& group, const QString& path, bool& updated) {
	int groupMax = 0; 
	auto vars = group.getVars();
	for (const auto& varPair : vars) {
		const NcVar& var = varPair.second;
		if (var.isNull() || var.getDimCount() < 1) {
			continue;
		}
		else if (var.getDimCount() <= 2) {
			size_t fc = var.getDim(0).getSize();
			QString name = path.isEmpty() ? QString::fromStdString(var.getName())
				: path + PATH_SEP + QString::fromStdString(var.getName());
			updated = updated || (_frameCounts[name] != static_cast<int>(fc));
			_frameCounts[name] = static_cast<int>(fc);
			groupMax = qMax(groupMax, static_cast<int>(fc));
		}
	}

	// Store per-group max frame count
	_groupMaxFrameCount[path] = groupMax;

	auto groups = group.getGroups();
	for (const auto& groupPair : groups) {
		const NcGroup& childGroup = groupPair.second;
		QString childPath = path.isEmpty() ? QString::fromStdString(groupPair.first)
			: path + PATH_SEP + QString::fromStdString(groupPair.first);
		traverseAndUpdateGroups(childGroup, childPath, updated);
	}
}

bool NetcdfSource::initFile() {
    try {
        _ncfile = std::make_unique<NcFile>(_filename.toStdString(), NcFile::read);
    } catch (const exceptions::NcException& e) {
        qDebug() << _filename << ": failed to open in initFile()" << endl;
        return false;
    }

    qDebug() << _filename << ": building field list" << endl;
    _fieldList.clear();
    _scalarList.clear();
    _matrixList.clear();
    _frameCounts.clear();
	_groupMaxFrameCount.clear();

    // Use the new private method for traversal
    traverseGroups(*_ncfile, "");

    // Get strings (global attributes)
    _strings.clear();
    auto atts = _ncfile->getAtts();
    for (const auto& attPair : atts) {
        const NcAtt& att = attPair.second;
        try {
            std::string value;
            att.getValues(value);
            _strings[QString::fromStdString(att.getName())] = QString::fromStdString(value);
        } catch (...) {
            // Ignore non-string attributes
        }
    }

    setUpdateType(Timer);

    qDebug() << "netcdf file initialized";
    return true;
}

Kst::Object::UpdateType NetcdfSource::internalDataSourceUpdate() {
	// netcdf-cxx4 does not require explicit sync for read-only files, but if needed:
	// _ncfile->sync();

	bool updated = false;
	traverseAndUpdateGroups(*_ncfile, "", updated);

	// Refresh the variable cache in case the file has changed
	_varMap.clear();
	traverseAndCacheVars(*_ncfile, "");

	return updated ? Object::Updated : Object::NoChange;
}

int NetcdfSource::readScalar(double* v, const QString& field)
{
    if (!_ncfile) return 0;
    try {
        NcVar var = getVarByFullPath(field);
        if (!var.isNull() && var.getDimCount() == 0) {
            var.getVar(v);
            return 1;
        }
    }
    catch (...) {}
    return 0;
}

int NetcdfSource::readString(QString* stringValue, const QString& stringName)
{
	if (!_ncfile) return 0;
	try {
		NcGroupAtt att = _ncfile->getAtt(stringName.toStdString());
		if (!att.isNull()) {
			std::string value;
			att.getValues(value);
			*stringValue = QString::fromStdString(value);
			return 1;
		}
	}
	catch (...) {}
	return 0;
}

int NetcdfSource::readField(double* v, const QString& field, int s, int n) {
    if (!_ncfile) return 0;
    if (isIndexField(field)) {
        int count = n < 0 ? 1 : n;
        for (int i = 0; i < count; ++i) {
            v[i] = double(s + i);
        }
        return count;
    }
    try {
        NcVar var = getVarByFullPath(field);
        if (var.isNull() || var.getDimCount() != 1) {
            return -1;
        }
        size_t dimSize = var.getDim(0).getSize();
        if (static_cast<size_t>(s) >= dimSize) {
            return 0;
        }

		NcType varType = var.getType();
		NcType::ncType t = varType.getTypeClass();
		if (isNumeric(t)) {
			size_t count = n < 0 ? 1 : qMin(dimSize-static_cast<size_t>(s), static_cast<size_t>(n));
			std::vector<size_t> start{ static_cast<size_t>(s) };
			std::vector<size_t> cnt{ count };

            // Read all requested values as doubles
            var.getVar(start, cnt, v);
			if (isNumeric(t) && iv->metaScalars(field).contains("add_offset") && iv->metaScalars(field).contains("scale_factor")) {
				double add_offset = 1.0, scale_factor = 1.0;
				add_offset = iv->metaScalars(field)["add_offset"];
				scale_factor = iv->metaScalars(field)["scale_factor"];
				for (size_t i = 0; i < count; ++i) {
					v[i] = v[i] * scale_factor + add_offset;
				}
			}
            return static_cast<int>(count);
        }
		return -1;
    }
    catch (...) {
        return -1;
    }
}

int NetcdfSource::readMatrix(double* v, const QString& field)
{
    if (!_ncfile) return -1;
    try {
        NcVar var = getVarByFullPath(field);
        if (var.isNull() || var.getDimCount() != 2) {
            return -1;
        }
        size_t xSize = var.getDim(0).getSize();
        size_t ySize = var.getDim(1).getSize();
        std::vector<size_t> start{ 0, 0 };
        std::vector<size_t> count{ xSize, ySize };
        var.getVar(start, count, v);
        return static_cast<int>(xSize * ySize);
    }
    catch (...) {
        return -1;
    }
}

int NetcdfSource::samplesPerFrame(const QString& field) {
    if (!_ncfile) return 0;
	if (isIndexField(field)) {
		return 1;
	}
	try {
        NcVar var = getVarByFullPath(field);
        if (var.isNull() || var.getDimCount() != 1) {
            return 0;
        }
        return static_cast<int>(1);
    }
    catch (...) {
        return 0;
    }
}

int NetcdfSource::frameCount(const QString& field) const {
	if (field.isEmpty() || isIndexField(field)) {
		QString groupPath = groupPathForIndexField(field);
		if (_groupMaxFrameCount.contains(groupPath)) {
			return _groupMaxFrameCount.value(groupPath, 0);
		}
		return 0; 
	}
	else {
		return _frameCounts.value(field, 0);
	}
}

QString NetcdfSource::fileType() const {
	return "netCDF";
}

bool NetcdfSource::isEmpty() const {
	return frameCount() < 1;
}

const QString& NetcdfSource::typeString() const
{
	return netCdfTypeString;
}

const QString NetcdfSource::netcdfTypeKey()
{
	return ::netCdfTypeString;
}

bool NetcdfSource::isIndexField(const QString& field) 
{
	return field.endsWith(INDEX_SUFFIX) || field == INDEX_FIELD;
}

QString NetcdfSource::groupPathForIndexField(const QString& field) 
{
	if (field == INDEX_FIELD) return "";
	if (field.endsWith(INDEX_SUFFIX)) return field.left(field.length() - 6); 
	return QString();
}

bool NetcdfSource::isNumeric(netCDF::NcType::ncType type)
{
	switch (type) {
	case NC_BYTE:
	case NC_UBYTE: 
	case NC_SHORT:
	case NC_USHORT: 
	case NC_INT:
	case NC_UINT: 
	case NC_INT64:
	case NC_UINT64: 
	case NC_FLOAT: 
	case NC_DOUBLE:
		return true;
		break;
	default:
		return false;
	}
}

bool NetcdfSource::isNumericInt(netCDF::NcType::ncType type)
{
	switch (type) {
	case NC_BYTE:
	case NC_UBYTE:
	case NC_SHORT:
	case NC_USHORT:
	case NC_INT:
	case NC_UINT:
	case NC_INT64:
	case NC_UINT64:
		return true;
		break;
	default:
		return false;
	}
}

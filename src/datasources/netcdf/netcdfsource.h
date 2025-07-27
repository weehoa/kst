/***************************************************************************
                netcdf_source.h  -  netCDF data source reader
                             -------------------
    begin                : 28/01/2005
    copyright            : (C) 2004 Nicolas Brisset <nicodev@users.sourceforge.net>
    email                : kst@kde.org
    modified             : 03/16/05 by K. Scott
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef NETCDFSOURCE_H
#define NETCDFSOURCE_H

#include "datasource.h"
#include "dataplugin.h"

#include <netcdf>
#include <memory>

class DataInterfaceNetCdfScalar;
class DataInterfaceNetCdfString;
class DataInterfaceNetCdfVector;
class DataInterfaceNetCdfMatrix;

class NetcdfSource : public Kst::DataSource {
  public:
    NetcdfSource(Kst::ObjectStore *store, QSettings *cfg, const QString& filename, const QString& type, const QDomElement &element);

    ~NetcdfSource();

    bool initFile();

    Kst::Object::UpdateType internalDataSourceUpdate();

    virtual const QString& typeString() const;

    static const QString netcdfTypeKey();

    int readScalar(double *v, const QString& field);

    int readString(QString *stringValue, const QString& stringName);

    int readField(double *v, const QString& field, int s, int n);

    int readMatrix(double *v, const QString& field);

    int samplesPerFrame(const QString& field);

    int frameCount(const QString& field = QString()) const;

    QString fileType() const;

    bool isEmpty() const;

    void reset();

  private:
    QMap<QString, int> _frameCounts;

    std::unique_ptr<netCDF::NcFile> _ncfile; // Updated to use netCDF4 C++ library

    QMap<QString, QString> _strings;

    QStringList _scalarList;
    QStringList _fieldList;
    QStringList _matrixList;

    QMap<QString, int> _groupMaxFrameCount;
    QMap<QString, netCDF::NcVar> _varMap;

    void traverseGroups(const netCDF::NcGroup& group, const QString& path);
    void traverseAndCacheVars(const netCDF::NcGroup& group, const QString& path);
    void traverseAndUpdateGroups(const netCDF::NcGroup& group, const QString& path, bool& updated);
    netCDF::NcVar getVarByFullPath(const QString& fullPath) const;

    static bool isIndexField(const QString& field);
    static QString groupPathForIndexField(const QString& field);
    static bool isNumeric(netCDF::NcType::ncType type);
    static bool isNumericInt(netCDF::NcType::ncType type);

    friend class DataInterfaceNetCdfScalar;
    friend class DataInterfaceNetCdfString;
    friend class DataInterfaceNetCdfVector;
    friend class DataInterfaceNetCdfMatrix;
    DataInterfaceNetCdfScalar* is;
    DataInterfaceNetCdfString* it;
    DataInterfaceNetCdfVector* iv;
    DataInterfaceNetCdfMatrix* im;
};

#endif


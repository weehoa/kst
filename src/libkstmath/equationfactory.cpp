/***************************************************************************
 *                                                                         *
 *   copyright : (C) 2007 The University of Toronto                        *
 *                   netterfield@astro.utoronto.ca                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "equationfactory.h"

#include "debug.h"
#include "equation.h"
#include "datacollection.h"
#include "objectstore.h"

namespace Kst {

EquationFactory::EquationFactory()
: ObjectFactory() {
  registerFactory(Equation::staticTypeTag, this);
}


EquationFactory::~EquationFactory() {
}


DataObjectPtr EquationFactory::generateObject(ObjectStore *store, QXmlStreamReader& xml) {
  QString tag, expression, xVector, output, descriptiveName;
  bool interpolate = false;
  while (!xml.atEnd()) {
    const QString n = xml.name().toString();
    if (xml.isStartElement()) {
      if (n == "equation") {
        QXmlStreamAttributes attrs = xml.attributes();
        tag = attrs.value("tag").toString();
        expression = attrs.value("expression").toString();
        interpolate = attrs.value("interpolate").toString().toLower() == "true";
        xVector = attrs.value("xvector").toString();
        output = attrs.value("output").toString();
        if (attrs.value("descriptiveNameIsManual").toString() == "true") {
          descriptiveName = attrs.value("descriptiveName").toString();
        }
        Object::processShortNameIndexAttributes(attrs);

        //FIXME Don't know if this is a bug in QXmlStreamReader or what, but readElementNext takes us
        //past the </equation> node to the </objects> node...
        //eq = xml.readElementText();
      } else {
        return 0;
      }
    } else if (xml.isEndElement()) {
      if (n == "equation") {
        break;
      } else {
        Debug::self()->log(QObject::tr("Error creating equation from Kst file."), Debug::Warning);
        return 0;
      }
    }
    xml.readNext();
  }

  if (xml.hasError()) {
    return 0;
  }

  VectorPtr vector = 0;
  if (store && !xVector.isEmpty()) {
    vector = kst_cast<Vector>(store->retrieveObject(xVector));
  }

  if (!vector) {
    Debug::self()->log(QObject::tr("Error creating equation from Kst file.  Could not find xVector."), Debug::Warning);
    return 0;
  }

  EquationPtr equation = store->createObject<Equation>();
  Q_ASSERT(equation);

  equation->setEquation(expression);
  equation->setExistingXVector(vector, interpolate);
  equation->setDescriptiveName(descriptiveName);

  equation->writeLock();
  equation->registerChange();
  equation->unlock();

  return static_cast<DataObjectPtr>(equation);
}

}

// vim: ts=2 sw=2 et

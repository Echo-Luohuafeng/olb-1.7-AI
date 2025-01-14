/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Mathias J. Krause, Benjamin Förster
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/


#ifndef SERIALIZER_HH
#define SERIALIZER_HH

#include <algorithm>
#include <iostream>
#include <ostream>
#include <fstream>
#include "serializer.h"
#include "communication/mpiManager.h"
#include "core/singleton.h"
#include "olbDebug.h"
#include "io/fileName.h"
#include "io/serializerIO.h"


namespace olb {

////////// class Serializer //////////////////

Serializer::Serializer(Serializable& serializable, std::string fileName)
  : _serializable(serializable), _iBlock(0), _size(0), _fileName(fileName)
{ }


void Serializer::resetCounter()
{
  _iBlock = 0;
}

std::size_t Serializer::getSize() const
{
  return _size;
}

bool* Serializer::getNextBlock(std::size_t& sizeBlock, bool loadingMode)
{
  return _serializable.getBlock(_iBlock++, sizeBlock, loadingMode);
}

template<bool includeLogOutputDir>
bool Serializer::load(std::string fileName, bool enforceUint)
{
  validateFileName(fileName);
  std::ifstream istr(getFullFileName<includeLogOutputDir>(fileName).c_str());
  if (istr) {
    istr2serializer(*this, istr, enforceUint);
    istr.close();
    _serializable.postLoad();
    return true;
  }
  else {
    return false;
  }
}

template<bool includeLogOutputDir>
bool Serializer::save(std::string fileName, bool enforceUint)
{
  validateFileName(fileName);

  // Determine binary size through `getSerializableSize()` method
  computeSize();

  std::ofstream ostr (getFullFileName<includeLogOutputDir>(fileName).c_str());
  if (ostr) {
    serializer2ostr(*this, ostr, enforceUint);
    ostr.close();
    return true;
  }
  else {
    return false;
  }
}

bool Serializer::load(const std::uint8_t* buffer)
{
  buffer2serializer(*this, buffer);
  _serializable.postLoad();
  return true;
}

bool Serializer::save(std::uint8_t* buffer)
{
  serializer2buffer(*this, buffer);
  return true;
}

void Serializer::computeSize(bool enforceRecompute)
{
  // compute size (only if it wasn't computed yet or is enforced)
  if (enforceRecompute || _size == 0) {
    _size = _serializable.getSerializableSize();
  }
}

void Serializer::validateFileName(std::string &fileName)
{
  if (fileName == "") {
    fileName = _fileName;
  }
  if (fileName == "") {
    fileName = "Serializable";
  }
}

template<bool includeLogOutputDir>
const std::string Serializer::getFullFileName(const std::string& fileName)
{
  if constexpr(includeLogOutputDir){
    return singleton::directories().getLogOutDir() + createParallelFileName(fileName) + ".dat";
  } else {
    return createParallelFileName(fileName) + ".dat";
  }
}



/////////////// Serializable //////////////////////////

template<bool includeLogOutputDir>
bool Serializable::save(std::string fileName, const bool enforceUint)
{
  Serializer tmpSerializer(*this, fileName);
  return tmpSerializer.save<includeLogOutputDir>();
}

template<bool includeLogOutputDir>
bool Serializable::load(std::string fileName, const bool enforceUint)
{
  Serializer tmpSerializer(*this, fileName);
  return tmpSerializer.load<includeLogOutputDir>();
}

bool Serializable::save(std::uint8_t* buffer)
{
  Serializer tmpSerializer(*this);
  return tmpSerializer.save(buffer);
}

bool Serializable::load(const std::uint8_t* buffer)
{
  Serializer tmpSerializer(*this);
  return tmpSerializer.load(buffer);
}

}  // namespace olb

#endif



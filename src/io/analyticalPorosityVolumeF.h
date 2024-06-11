/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerländer, Dennis Teutscher
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

#ifndef ANALYTICAL_POROSITY_VOLUME_F_H
#define ANALYTICAL_POROSITY_VOLUME_F_H

#ifdef FEATURE_VDB
#include <openvdb/openvdb.h>
#include <openvdb/tree/ValueAccessor.h>
#include <openvdb/Grid.h>
#endif

#if defined(FEATURE_VTK)
#include <vtkNew.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#endif

namespace olb {

template <typename T>
class AnalyticalPorosityVolumeF final : public AnalyticalF<3,T,T> {
private:
#if defined(FEATURE_VTK)
  vtkNew<vtkStructuredPointsReader> _reader;
  vtkStructuredPoints* _data;
#endif
#ifdef FEATURE_VDB
  openvdb::GridPtrVecPtr grids;
#endif
  std::string _fileName;
  Vector<int,3> _shape;
  T _spacing;

public:
  AnalyticalPorosityVolumeF(std::string fileName, T spacing)
      : AnalyticalF<3,T,T>(1), _spacing{spacing}
  {
    OstreamManager clout(std::cout, "PorosityVolumeImporter");
    _fileName = fileName;
    #ifdef FEATURE_VDB
    if (isVDBFile(_fileName)) {
      openvdb::initialize();
      openvdb::io::File _vdbFile(_fileName);
      // Read the OpenVDB file
      _vdbFile.open();
      grids = _vdbFile.getGrids();
      auto grid = openvdb::gridPtrCast<openvdb::FloatGrid>(*grids->begin());
      _vdbFile.close();
      openvdb::CoordBBox bounds = grid->evalActiveVoxelBoundingBox();
      _shape[0] = bounds.max().x() - bounds.min().x();
      _shape[1] = bounds.max().y() - bounds.min().y();
      _shape[2] = bounds.max().z() - bounds.min().z();
    }
    #endif
    #if not defined(FEATURE_VTK)
    if (isVTKFile(_fileName))
    {
      std::cerr << "To use the VTK format, add VTK to the FEATURES list in config.mk";
      exit(1);
    }
    #endif
    #if defined(FEATURE_VTK)
    if (isVTKFile(_fileName)) {
      _reader->SetFileName(fileName.c_str());
      _reader->Update();
      _reader->SetScalarsName(_reader->GetScalarsNameInFile(0));
      auto* _data = _reader->GetOutput();
      _data->GetDimensions(_shape.data());
    }
    #endif
    #ifndef FEATURE_VDB
    if (isVDBFile(_fileName))
    {
      std::cerr << "To use the VDB format, add VDB to the FEATURES list in config.mk";
      exit(1);
    }
    #endif
  }

  Vector<T,3> getPhysShape() const {
    return _shape * _spacing;
  }

  bool operator()(T output[], const T physR[]) override{
    #ifdef FEATURE_VDB
    if (isVDBFile(_fileName)) {
      int iX;
      int iY;
      int iZ;
      iX = util::floor(physR[0] / _spacing);
      iY = util::floor(physR[1] / _spacing);
      iZ = util::floor(physR[2] / _spacing);
      if (iX >= 0 && iY >= 0 && iZ >= 0 && iX < _shape[0] && iY < _shape[1] && iZ < _shape[2]) {
        auto grid = openvdb::gridPtrCast<openvdb::FloatGrid>(*grids->begin());
        openvdb::Coord location(iX, iY, iZ);
        openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
        output[0] =accessor.getValue(location);
      } else {
        output[0]=0;
      }
      return true;
    }
    #endif
    #if defined(FEATURE_VTK)
    if(isVTKFile(_fileName)) {
      int iX;
      int iY;
      int iZ;
      iX = util::floor(physR[0] / _spacing);
      iY = util::floor(physR[1] / _spacing);
      iZ = util::floor(physR[2] / _spacing);
      if (iX >= 0 && iY >= 0 && iZ >= 0 && iX < _shape[0] && iY < _shape[1] && iZ < _shape[2]) {
        auto* _data = _reader->GetOutput();
        output[0] = *static_cast<T*>(_data->GetScalarPointer(iX, iY, iZ));
      } else {
        output[0] = 0;
      }
      return true;
    }
    #endif
    return false;
  }

private:
  bool isVDBFile(const std::string& fileName) {
    return (fileName.size() >= 4 && fileName.substr(fileName.size() - 4) == ".vdb");
  }
  bool isVTKFile(const std::string& fileName) {
    return (fileName.size() >= 4 && fileName.substr(fileName.size() - 4) == ".vtk");
  }

};

}

#endif

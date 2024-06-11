/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Albert Mink, Mathias J. Krause
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

/** \file
 * A method to write vtk data for cuboid geometries
 * (only for uniform grids) -- generic implementation.
 */

#ifndef SUPER_VTM_WRITER_2D_HH
#define SUPER_VTM_WRITER_2D_HH

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "core/singleton.h"
#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry2D.h"
#include "communication/mpiManager.h"
#include "io/base64.h"
#include "io/fileName.h"
#include "io/superVtmWriter2D.h"


namespace olb {


template<typename T, typename OUT_T, typename W>
SuperVTMwriter2D<T,OUT_T,W>::SuperVTMwriter2D( std::string name, int overlap, bool binary )
  : clout( std::cout,"SuperVTMwriter2D" ), _createFile(false), _name(name), _overlap(overlap), _binary(binary)
{
  static_assert(std::is_same_v<OUT_T, float> || std::is_same_v<OUT_T, double>,
              "OUT_T must be either float or double");
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::write(int iT)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  //  !!!!!!!!!!! check whether _pointerVec is empty
  if ( _pointerVec.empty() ) {
    clout << "Error: Did you add a Functor ?";
  }
  else {
    // no gaps between vti files (cuboids)
    for (SuperF2D<T,W>* f : _pointerVec) {
      f->getSuperStructure().communicate();
    }

    // to get first element _pointerVec
    // problem if functors with different SuperStructure are stored
    // since till now, there is only one origin
    auto it_begin = _pointerVec.cbegin();
    if (it_begin == _pointerVec.end()) {
      throw std::runtime_error("No functor to write");
    }
    CuboidGeometry2D<T> const& cGeometry = (**it_begin).getSuperStructure().getCuboidGeometry();
    // no gaps between vti files (cuboids)
    LoadBalancer<T>& load = (**it_begin).getSuperStructure().getLoadBalancer();

    // PVD, owns all
    if ( rank == 0 ) {
      std::string pathPVD = singleton::directories().getVtkOutDir()
                            + createFileName( _name ) + ".pvd";
      dataPVDmaster( iT, pathPVD,  "data/" + createFileName( _name, iT ) + ".vtm" );

      std::string pathVTM = singleton::directories().getVtkOutDir()
                            + "data/" + createFileName( _name, iT ) + ".vtm";
      preambleVTM(pathVTM);
      for (int iC = 0; iC < cGeometry.getNc(); iC++) {
        dataVTM( iC, pathVTM, createFileName( _name, iT, iC) + ".vti" );
      }
      closeVTM(pathVTM);
    }
    // VTI, each process writes his cuboids
    int originLatticeR[3] = {int()};
    for (int iCloc = 0; iCloc < load.size(); iCloc++) {
      int nx = cGeometry.get(load.glob(iCloc)).getNx();
      int ny = cGeometry.get(load.glob(iCloc)).getNy();
      // to be changed into the following line once local refinement has been implemented
      // double deltaX = cGeometry.get(load.glob(iCloc)).getDeltaR();
      T delta = cGeometry.getMotherCuboid().getDeltaR();

      std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                                + createFileName( _name, iT, load.glob(iCloc) ) + ".vti";

      // get dimension/extent for each cuboid
      originLatticeR[0] = load.glob(iCloc);
      T originPhysR[2] = {T()};
      cGeometry.getPhysR(originPhysR,originLatticeR);

      preambleVTI(fullNameVTI, -_overlap,-_overlap, nx+_overlap-1, ny+_overlap-1, originPhysR[0],originPhysR[1], delta);
      for (auto it : _pointerVec) {
        if (_binary) {
          dataArrayBinary(fullNameVTI, (*it), load.glob(iCloc), nx,ny);
        }
        else {
          dataArray(fullNameVTI, (*it), load.glob(iCloc), nx,ny);
        }
      }
      closePiece(fullNameVTI);
      closeVTI(fullNameVTI);
    }
  }
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T, W>::write(SuperF2D<T,W>& f, int iT)
{
  CuboidGeometry2D<T> const& cGeometry = f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();
  // no gaps between vti files (cuboids)
  f.getSuperStructure().communicate();
  T delta = cGeometry.getMotherCuboid().getDeltaR();

  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  // write a pvd file, which links all vti files
  // each vti file is written by one thread, which may own severals cuboids
  if ( rank == 0 ) {
    // master only
    std::string pathVTM = singleton::directories().getVtkOutDir()
                          + createFileName( f.getName(), iT )  + ".vtm";

    preambleVTM(pathVTM);
    for (int iC = 0; iC < cGeometry.getNc(); iC++) {
      std::string nameVTI = "data/" + createFileName( f.getName(), iT, iC) + ".vti";
      // puts name of .vti piece to a .pvd file [fullNamePVD]
      dataVTM( iC, pathVTM, nameVTI );
    }
    closeVTM(pathVTM);
  } // master only

  for (int iCloc = 0; iCloc < load.size(); iCloc++) {
    // cuboid
    int nx = cGeometry.get(load.glob(iCloc)).getNx();
    int ny = cGeometry.get(load.glob(iCloc)).getNy();
    // to be changed into the following line once local refinement has been implemented
    // double deltaX = cGeometry.get(load.glob(iCloc)).getDeltaR();

    std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName( f.getName(), iT, load.glob(iCloc) ) + ".vti";

    // get dimension/extent for each cuboid
    int const originLatticeR[3] = {load.glob(iCloc),0,0};
    T originPhysR[2] = {T()};
    cGeometry.getPhysR(originPhysR,originLatticeR);

    preambleVTI(fullNameVTI, -_overlap,-_overlap, nx+_overlap-1, ny+_overlap-1, originPhysR[0],originPhysR[1], delta);
    if (_binary) {
      dataArrayBinary(fullNameVTI, f, load.glob(iCloc), nx,ny);
    }
    else {
      dataArray(fullNameVTI, f, load.glob(iCloc), nx,ny);
    }
    closePiece(fullNameVTI);
    closeVTI(fullNameVTI);
  } // cuboid
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::write(std::shared_ptr<SuperF2D<T,W>> ptr_f, int iT)
{
  write(*ptr_f, iT);
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::createMasterFile()
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if ( rank == 0 ) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName( _name ) + ".pvd";
    preamblePVD(fullNamePVDmaster);
    closePVD(fullNamePVDmaster);
    _createFile = true;
  }
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::addFunctor(SuperF2D<T,W>& f)
{
  _pointerVec.push_back(&f);
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::addFunctor(SuperF2D<T,W>& f, const std::string& functorName)
{
  f.getName() = functorName;
  _pointerVec.push_back(&f);
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::clearAddedFunctors()
{
  _pointerVec.clear();
}

template<typename T, typename OUT_T, typename W>
std::string SuperVTMwriter2D<T,OUT_T,W>::getName() const
{
  return _name;
}




////////////////////private member functions///////////////////////////////////
template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::preambleVTI (const std::string& fullName,
    int x0, int y0, int x1, int y1, T originX, T originY, T delta)
{
  const BaseType<T> d_delta = delta;
  const BaseType<T> d_origin[2] = {originX, originY};
  std::ofstream fout(fullName, std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n";
  fout << "<ImageData WholeExtent=\""
       << x0 <<" "<< x1 <<" "
       << y0 <<" "<< y1 <<" "
       << 0 <<" "<< 0
       << "\" Origin=\"" << d_origin[0] << " " << d_origin[1] << " " << "0"
       << "\" Spacing=\"" << d_delta << " " << d_delta << " " << d_delta << "\">\n";
  fout << "<Piece Extent=\""
       << x0 <<" "<< x1 <<" "
       << y0 <<" "<< y1 <<" "
       << 0  <<" "<< 0 <<"\">\n";
  fout << "<PointData>\n";
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::closeVTI(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece, std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</ImageData>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::preamblePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD, std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n"
       << "<Collection>\n";
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::closePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD, std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::preambleVTM(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD, std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
       << "byte_order=\"LittleEndian\">\n"
       << "<vtkMultiBlockDataSet>\n" ;
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::closeVTM(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD, std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</vtkMultiBlockDataSet>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::dataVTM(int iC, const std::string& fullNamePVD,
                                    const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD, std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<Block index=\"" << iC << "\" >\n";
  fout << "<DataSet index= \"0\" " << "file=\"" << namePiece << "\">\n"
       << "</DataSet>\n";
  fout << "</Block>\n";
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::dataPVDmaster(int iT,
    const std::string& fullNamePVDMaster,
    const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster, std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25,std::ios::end);    // jump -25 form the end of file to overwrite closePVD

    fout << "<DataSet timestep=\"" << iT << "\" "
         << "group=\"\" part=\"\" "
         << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  }
  else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::dataArray(const std::string& fullName,
                                      SuperF2D<T,W>& f, int iC, int nx, int ny)
{
  std::cout << "DIOCANE" <<std::endl <<std::endl <<std::endl <<std::endl <<std::endl;
  std::ofstream fout( fullName, std::ios::out | std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  fout << "<DataArray " ;
  if constexpr (std::is_same_v<OUT_T, float>) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "NumberOfComponents=\"" << f.getTargetDim() <<"\">\n";
  }
  else if constexpr (std::is_same_v<OUT_T, double>) {
    fout << "type=\"Float64\" Name=\"" << f.getName() << "\" "
         << "NumberOfComponents=\"" << f.getTargetDim() <<"\">\n";
  }

  int i[3] = {iC, 0, 0};
  W evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = W();
  }
  // since cuboid has been blowed up by _overlap [every dimension]
  // looping from -_overlap to ny (= ny+_overlap, as passed)
  std::vector<int> tmpVec( 3,int(0) );
  for (i[2]=-_overlap; i[2] < ny+_overlap; ++i[2]) {
    for (i[1]=-_overlap; i[1] < nx+_overlap; ++i[1]) {
      f(evaluated,i);
      for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
        //  tmpVec = {iC,iX,iY,iZ};  // std=c++11
        fout <<  evaluated[iDim] << " ";
      }
    }
  }
  fout << "</DataArray>\n";
  fout.close();
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::dataArrayBinary(const std::string& fullName,
    SuperF2D<T,W>& f, int iC, int nx, int ny)
{
  std::ofstream fout( fullName, std::ios::out | std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  fout << "<DataArray ";
  if constexpr (std::is_same_v<OUT_T, float>) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "format=\"binary\" encoding=\"base64\" "
         << "NumberOfComponents=\"" << f.getTargetDim() <<"\">\n";
  }
  else if constexpr (std::is_same_v<OUT_T, double>) {
    fout << "type=\"Float64\" Name=\"" << f.getName() << "\" "
         << "format=\"binary\" encoding=\"base64\" "
         << "NumberOfComponents=\"" << f.getTargetDim() <<"\">\n";
  }
  fout.close();

  std::ofstream ofstr( fullName, std::ios::out | std::ios::app | std::ios::binary );
  if (!ofstr) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  size_t fullSize = f.getTargetDim() * (nx+2*_overlap) * (ny+2*_overlap);    //  how many numbers to write
  size_t binarySize = size_t( fullSize*sizeof(OUT_T) );
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int uintBinarySize = (unsigned int)binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<OUT_T>* dataEncoder = new Base64Encoder<OUT_T>( ofstr, fullSize );

  int i[3] = {iC, 0, 0};
  W evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = W();
  }
  int itter = 0;
  std::unique_ptr<OUT_T[]> bufferFloat(new OUT_T[fullSize]);
  for (i[2] = -_overlap; i[2] < ny+_overlap; ++i[2]) {
    for (i[1] = -_overlap; i[1] < nx+_overlap; ++i[1]) {
      f(evaluated,i);
      for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
        bufferFloat[ itter ] = OUT_T( evaluated[iDim] );
        itter++;
      }
    }
  }
  dataEncoder->encode( &bufferFloat[0], fullSize );
  ofstr.close();

  std::ofstream ffout( fullName,  std::ios::out | std::ios::app );
  if (!ffout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  ffout << "\n</DataArray>\n";
  ffout.close();
  delete dataEncoder;
}

template<typename T, typename OUT_T, typename W>
void SuperVTMwriter2D<T,OUT_T,W>::closePiece(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece, std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</PointData>\n";
  fout << "</Piece>\n";
  fout.close();
}


}  // namespace olb

#endif

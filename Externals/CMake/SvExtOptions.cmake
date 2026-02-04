# Copyright (c) 2014-2015 The Regents of the University of California.
# All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#-----------------------------------------------------------------------------
# Toplevel directories for src, bin, build or externals
set(SV_EXTERNALS_TOPLEVEL_DIR "${CMAKE_BINARY_DIR}/svExternals")

message(STATUS "[SvExtOptions] ")
message(STATUS "[SvExtOptions] --------- SvExtOptions ----------")

set(SV_EXTERNALS_TOPLEVEL_SRC_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/src"
  CACHE PATH "Directory where source files for externals will be put")
set(SV_EXTERNALS_TOPLEVEL_BIN_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/bin/${SV_COMPILER_DIR}/${SV_COMPILER_VERSION_DIR}/${SV_ARCH_DIR}/${SV_BUILD_TYPE_DIR}"
  CACHE PATH "Directory where install files for externals will be put")
set(SV_EXTERNALS_TOPLEVEL_BLD_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/build"
  CACHE PATH "Directory where build files for externals will be put")
set(SV_EXTERNALS_TOPLEVEL_PFX_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/prefix"
  CACHE PATH "Directory where prefix files for externals will be put")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Tar install directory
simvascular_today(YEAR MONTH DAY)
option(SV_EXTERNALS_TAR_BINARIES "Tar up the final built or downloaded externals" OFF)
set(SV_EXTERNALS_TAR_INSTALL_DIR "${CMAKE_BINARY_DIR}/tar_output/${YEAR}.${MONTH}.${DAY}"
  CACHE PATH "Directory where source files for externals will be put")
file(MAKE_DIRECTORY "${SV_EXTERNALS_TAR_INSTALL_DIR}")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# URLs for external downloads and git repositories
set(SV_EXTERNALS_URL "" CACHE STRING "SimVascular Externals base URL (set this to a mirror)")
set(SV_EXTERNALS_ORIGINALS_URL "" CACHE STRING "URL with source downloads for externals")
if(SV_EXTERNALS_URL AND NOT SV_EXTERNALS_ORIGINALS_URL)
  set(SV_EXTERNALS_ORIGINALS_URL "${SV_EXTERNALS_URL}/${SV_EXTERNALS_VERSION_NUMBER}/src/originals" CACHE STRING "URL with source downloads for externals" FORCE)
endif()
if(NOT SV_EXTERNALS_URL AND NOT SV_EXTERNALS_ORIGINALS_URL)
  message(WARNING "SV_EXTERNALS_URL is unset. Set it to a valid mirror, set SV_EXTERNALS_ORIGINALS_URL directly, or provide per-external SV_EXTERNALS_<proj>_MANUAL_SOURCE_URL values.")
endif()

if(SV_EXTERNALS_VERSION_NUMBER VERSION_EQUAL "2024.05" AND NOT SV_EXTERNALS_URL AND NOT SV_EXTERNALS_ORIGINALS_URL)
  if(NOT SV_EXTERNALS_QT_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_QT_MANUAL_SOURCE_URL "https://download.qt.io/official_releases/qt/6.6/6.6.2/single/qt-everywhere-src-6.6.2.tar.xz" CACHE STRING "Manual specification of QT source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_HDF5_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_HDF5_MANUAL_SOURCE_URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.gz" CACHE STRING "Manual specification of HDF5 source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_TINYXML2_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_TINYXML2_MANUAL_SOURCE_URL "https://github.com/leethomason/tinyxml2/archive/refs/tags/6.2.0.tar.gz" CACHE STRING "Manual specification of TINYXML2 source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_PYTHON_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_PYTHON_MANUAL_SOURCE_URL "https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tgz" CACHE STRING "Manual specification of PYTHON source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_PIP_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_PIP_MANUAL_SOURCE_URL "https://bootstrap.pypa.io/get-pip.py" CACHE STRING "Manual specification of PIP source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_NUMPY_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_NUMPY_MANUAL_SOURCE_URL "https://github.com/numpy/numpy/releases/download/v1.14.3/numpy-1.14.3.tar.gz" CACHE STRING "Manual specification of NUMPY source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_FREETYPE_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_FREETYPE_MANUAL_SOURCE_URL "https://download.savannah.gnu.org/releases/freetype/freetype-2.13.0.tar.gz" CACHE STRING "Manual specification of FREETYPE source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_SWIG_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_SWIG_MANUAL_SOURCE_URL "https://downloads.sourceforge.net/swig/swig-3.0.12.tar.gz" CACHE STRING "Manual specification of SWIG source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_MMG_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_MMG_MANUAL_SOURCE_URL "https://github.com/MmgTools/mmg/archive/refs/tags/v5.3.9.tar.gz" CACHE STRING "Manual specification of MMG source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_GDCM_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_GDCM_MANUAL_SOURCE_URL "https://github.com/malaterre/GDCM/archive/refs/tags/v3.0.10.tar.gz" CACHE STRING "Manual specification of GDCM source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_VTK_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_VTK_MANUAL_SOURCE_URL "https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz" CACHE STRING "Manual specification of VTK source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_ITK_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_ITK_MANUAL_SOURCE_URL "https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.4.0/InsightToolkit-5.4.0.tar.gz" CACHE STRING "Manual specification of ITK source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_OpenCASCADE_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_OpenCASCADE_MANUAL_SOURCE_URL "https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V7_6_0.tar.gz" CACHE STRING "Manual specification of OpenCASCADE source URL" FORCE)
  endif()
  if(NOT SV_EXTERNALS_MITK_MANUAL_SOURCE_URL)
    set(SV_EXTERNALS_MITK_MANUAL_SOURCE_URL "https://github.com/MITK/MITK/archive/refs/tags/v2022.10.tar.gz" CACHE STRING "Manual specification of MITK source URL" FORCE)
  endif()
endif()
mark_as_advanced(SV_EXTERNALS_ORIGINALS_URL)
set(SV_EXTERNALS_BINARIES_URL_PREFIX "${SV_PLATFORM_DIR}.${SV_COMPILER_DIR}-${SV_COMPILER_VERSION_DIR}.${SV_ARCH_DIR}" CACHE STRING "String that gets appended to the beginning of each individual external pre-built binary download")
mark_as_advanced(SV_EXTERNALS_BINARIES_URL_PREFIX)
set(SV_EXTERNALS_GIT_URL "http://github.com/SimVascular" CACHE STRING "Git URL for SimVascular")
mark_as_advanced(SV_EXTERNALS_GIT_URL)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# SV_EXTERNALS_MAC_PACKAGE_MANAGER
if (APPLE)
  message(WARNING "Need to make sure to have openssl from homebrew or macports. Specify which is used with SV_EXTERNALS_MAC_PACKAGE_MANAGER variable.")
  set(SV_EXTERNALS_MAC_PACKAGE_MANAGER "HOMEBREW" CACHE STRING "Options are HOMEBREW OR MACPORTS")
  set_property(CACHE SV_EXTERNALS_MAC_PACKAGE_MANAGER PROPERTY STRINGS HOMEBREW MACPORTS)
endif()
#-----------------------------------------------------------------------------

# Add externals with default values of version, build_with, shared, dirname,
# and optional install dirname. Order matters; put independent packages first
# Must have existing "EXTERNAL_NAME.cmake" file underneath CMake
# (i.e. QT.cmake for QT)
# "EXTERNAL_NAME" "ENABLE_EXTERNAL" "BUILD_SHARED" "BUILD_DIR_NAME" "INSTALL_DIR_NAME"
#-----------------------------------------------------------------------------
# QT
message(STATUS "[SvExtOptions] QT sv_externals_add_new_external ")
#option(SV_QT_DIR "The location of the Qt package install directory" "")
#option(SV_EXTERNALS_USE_PREBUILT_QT "Instead of downloading or building, use a specified QT" OFF)

sv_externals_add_new_external(QT ${SV_EXTERNALS_QT_VERSION} ON ON qt qt)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ML
if (SV_EXTERNALS_VERSION_NUMBER GREATER_EQUAL "2019.02")
  sv_externals_add_new_external(ML ${SV_EXTERNALS_ML_VERSION} ON ON ml ml)
endif()

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# HDF5
if (SV_EXTERNALS_VERSION_NUMBER GREATER_EQUAL "2018.05")
  sv_externals_add_new_external(HDF5 ${SV_EXTERNALS_HDF5_VERSION} ON ON hdf5 hdf5)
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # TINYXML2
  sv_externals_add_new_external(TINYXML2 ${SV_EXTERNALS_TINYXML2_VERSION} ON ON tinyxml2 tinyxml2)
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# PYTHON
sv_externals_add_new_external(PYTHON ${SV_EXTERNALS_PYTHON_VERSION} ON ON python python)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# PIP
sv_externals_add_new_external(PIP ${SV_EXTERNALS_PIP_VERSION} ON ON pip none)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# NUMPY
sv_externals_add_new_external(NUMPY ${SV_EXTERNALS_NUMPY_VERSION} ON ON numpy none)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# FREETYPE
sv_externals_add_new_external(FREETYPE ${SV_EXTERNALS_FREETYPE_VERSION} ON ON freetype freetype)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# SWIG
sv_externals_add_new_external(SWIG ${SV_EXTERNALS_SWIG_VERSION} OFF ON swig swig)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# MMG
sv_externals_add_new_external(MMG ${SV_EXTERNALS_MMG_VERSION} ON OFF mmg mmg)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# GDCM
sv_externals_add_new_external(GDCM ${SV_EXTERNALS_GDCM_VERSION} ON ON gdcm gdcm)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# VTK
sv_externals_add_new_external(VTK ${SV_EXTERNALS_VTK_VERSION} ON ON vtk vtk)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ITK
sv_externals_add_new_external(ITK ${SV_EXTERNALS_ITK_VERSION} ON ON itk itk)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# OpenCASCADE
sv_externals_add_new_external(OpenCASCADE ${SV_EXTERNALS_OpenCASCADE_VERSION} ON ON opencascade opencascade)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# MITK
option(SV_EXTERNALS_BUILD_MITK_WITH_PYTHON "Build MITK with python" ON)
sv_externals_add_new_external(MITK ${SV_EXTERNALS_MITK_VERSION} ON ON mitk mitk)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#Download options for tcltk
option(SV_EXTERNALS_DOWNLOAD_TCLTK "Download instead of build TCLTK" ON)

#message(FATAL_ERROR "[SvExtOptions] SV_EXTERNALS_LIST: ${SV_EXTERNALS_LIST}")

message(STATUS "[SvExtOptions] ----- Done SvExtOptions -----")
message(STATUS "[SvExtOptions] ")


 # Tarang-4.0
 # 
 # Copyright (C) 2008, 2009  Mahendra K. Verma
 #
 # Mahendra K. Verma
 # Indian Institute of Technology, Kanpur-208016
 # UP, India
 #
 # mkv@iitk.ac.in
 #
 # This file is part of Tarang-4.0 .
 #
 # Tarang-4.0 is free software; you can redistribute it and/or
 # modify it under the terms of the GNU General Public License
 # as published by the Free Software Foundation; either version 2
 # of the License, or (at your option) any later version.
 # Tarang-4.0 is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU General Public License
 # along with Tarang-4.0; if not, write to the Free Software
 # Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 #
PROJECT(TARANG_MPI)
cmake_minimum_required(VERSION 2.6)

SUBDIRS(lib src)

ADD_DEFINITIONS(-DTRANSPOSE)

INCLUDE_DIRECTORIES(
${CMAKE_SOURCE_DIR}/lib/basis_basicfn
${CMAKE_SOURCE_DIR}/lib/scft 
${CMAKE_SOURCE_DIR}/lib/four 
${CMAKE_SOURCE_DIR}/lib/plainfields
${CMAKE_SOURCE_DIR}/lib/fields
${CMAKE_SOURCE_DIR}/lib/IncFlow 
${CMAKE_SOURCE_DIR}/lib/IncFluid
#/Users/supriyo/local/include
)


LINK_DIRECTORIES(
${CMAKE_SOURCE_DIR}/lib/basis_basicfn
${CMAKE_SOURCE_DIR}/lib/scft 
${CMAKE_SOURCE_DIR}/lib/four 
${CMAKE_SOURCE_DIR}/lib/plainfields
${CMAKE_SOURCE_DIR}/lib/fields
${CMAKE_SOURCE_DIR}/lib/IncFlow 
${CMAKE_SOURCE_DIR}/lib/IncFluid
${CMAKE_SOURCE_DIR}/src
${CMAKE_SOURCE_DIR}/src/IFluid
${CMAKE_SOURCE_DIR}/src/IMHD
${CMAKE_SOURCE_DIR}/src/IScalar
${CMAKE_SOURCE_DIR}/src/RBslip
#/Users/supriyo/local/lib
)


ADD_EXECUTABLE(tarangMPI  src/main.h src/main.cc  
	src/read_prog_para.cc
	src/read_para.cc 
	src/read_diag_para.cc
	src/IFluid/Ifluid_main.h  
	src/IFluid/Ifluid_main.cc 
	src/IFluid/Ifluid_diag_main.cc 
	src/IScalar/Iscalar_main.h  
	src/IScalar/Iscalar_main.cc
	src/IScalar/Iscalar_diag_main.cc
	src/IMHD/IMHD_main.h  
	src/IMHD/IMHD_main.cc
	src/IMHD/IMHD_diag_main.cc
	src/RBslip/RB_slip_main.h   
	src/RBslip/RB_slip_main.cc 
	src/RBslip/RB_slip_diag_main.cc
	)



TARGET_LINK_LIBRARIES(tarangMPI IncFluid IncFlow fields plainfields four scft basis_basic m blitz fftw3_mpi fftw3 mpi)

#find_package(HDF5 REQUIRED)
SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -g")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")


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
#PROJECT(tarangMPI)

SET(CMAKE_CXX_COMPILER mpicxx)

#ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK -Wall)

ADD_DEFINITIONS(-Wall)

#Choose compiler
#CC=gcc-3.3 CXX=g++-3.3 cmake

# Make sure the compiler can find include files from our lib library.
INCLUDE_DIRECTORIES(
${CMAKE_SOURCE_DIR}/lib/basis_basicfn
${CMAKE_SOURCE_DIR}/lib/scft 
${CMAKE_SOURCE_DIR}/lib/four 
${CMAKE_SOURCE_DIR}/lib/plainfields
${CMAKE_SOURCE_DIR}/lib/fields
${CMAKE_SOURCE_DIR}/lib/IncFlow 
${CMAKE_SOURCE_DIR}/lib/IncFluid
)

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(
${CMAKE_SOURCE_DIR}/lib/basis_basicfn
${CMAKE_SOURCE_DIR}/lib/scft 
${CMAKE_SOURCE_DIR}/lib/four 
${CMAKE_SOURCE_DIR}/lib/plainfields
${CMAKE_SOURCE_DIR}/lib/fields
${CMAKE_SOURCE_DIR}/lib/IncFlow 
${CMAKE_SOURCE_DIR}/lib/IncFluid
)

SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3 -pipe")

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx".  The extensions are automatically found.
ADD_EXECUTABLE(../tarangMPI  main.h main.cc  
	read_prog_para.cc
	read_para.cc 
	IFluid/Ifluid_main.h  
	IFluid/Ifluid_main.cc 
	IPassiveScalar/Ipassive_scalar_main.h  
	IPassiveScalar/Ipassive_scalar_main.cc
	IMHD/IMHD_main.h  
	IMHD/IMHD_main.cc
	RBslip/RB_slip_main.h   
	RBslip/RB_slip_main.cc   
	RBslip/read_RB_slip_para.cc)


# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(../tarangMPI IncFluid IncFlow fields plainfields four scft basis_basic m blitz fftw3_mpi fftw3 mpi)

SET(CMAKE_BUILD_TYPE distribution)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")

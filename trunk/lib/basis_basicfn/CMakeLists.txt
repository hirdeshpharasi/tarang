
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

# ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)

#INCLUDE_DIRECTORIES(/Users/mkv/local/include )

#LINK_DIRECTORIES(/Users/mkv/local/lib )

ADD_LIBRARY(basis_basic 
		blitz_random.h 
		extern_def_vars.h
		basis_basicfn_inline.h basis_basicfn.h 
		basis_basicfn.cc)

# TARGET_LINK_LIBRARIES(basis_basic fftw3_mpi fftw3  mpi m)

#SET(CMAKE_BUILD_TYPE distribution)
#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_CXX_FLAGS "-Wall -Werror")


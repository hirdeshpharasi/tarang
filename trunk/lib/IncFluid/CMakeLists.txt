
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
#SET(CMAKE_CXX_COMPILER mpicxx)

# ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)

# ADD_DEFINITIONS(-DWAVENOACTUAL -DANISDIRN1)

# INCLUDE_DIRECTORIES(/Users/mkv/local/include )

# LINK_DIRECTORIES(/Users/mkv/local/lib )

ADD_LIBRARY(IncFluid  IncFluid_constants.h 
	IncFluid_Time.h	IncFluid_Time.cc
	IncFluid.h IncFluid.cc 
	add_force.cc
	add_pressure_grad.cc 
	compute_dt.cc
	compute_rhs.cc 
	file_operations.cc
	force/compute_force_ABC.cc
	force/compute_force_const_ek_hk.cc
	force/compute_force_decay.cc
	force/compute_force_dynamo_6mode.cc
	force/compute_force_ek_hk_supply.cc
	force/compute_force_LM.cc
	force/compute_force_main.cc
	force/compute_force_modes.cc
	force/compute_force_RB.cc
	force/compute_force_TG.cc
	input/init_cond_ABC.cc
	input/init_cond_dynamo_6mode.cc
	input/init_cond_energy.cc
	input/init_cond_field.cc
	input/init_cond_main.cc
	input/init_cond_modes.cc
	input/init_cond_RB_misc.cc
	input/init_cond_TG.cc
	output/output_energy.cc
	output/output_ET.cc
	output/output_field.cc
	output/output_main.cc
	output/output_moment.cc
	output/output_str_fn.cc
	time_advance/single_time_step.cc
	time_advance/time_advance_RBslip.cc
	time_advance/time_advance_RBslipMHD.cc
    	time_advance/time_advance_SF.cc 
	time_advance/time_advance_VF_MHD.cc 
  	time_advance/time_advance.cc 
	)

# TARGET_LINK_LIBRARIES(IncFluid fftw3_mpi fftw3  mpi m)

#SET(CMAKE_BUILD_TYPE distribution)
#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_CXX_FLAGS "-Wall -Werror")

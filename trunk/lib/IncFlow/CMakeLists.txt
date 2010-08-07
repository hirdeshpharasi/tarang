
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

ADD_LIBRARY(IncFlow 
	IncVF 
	IncSF 
	Nlin 
	EnergyTr.h   
	compute_divergence.cc
	compute_pressure.cc
	force_spectrum.cc
	nlin_calc/RSprod.cc 
	nlin_calc/compute_derivative.cc	
	nlin_calc/compute_nlin.cc 
	nlin_calc/compute_nlin_SF.cc 
	nlin_calc/compute_nlin_VF_MHD.cc
	nlin_calc/compute_nlin_VF_SF.cc
	misc/add_nlin_dt.cc
	misc/copy_field.cc
	misc/dealias.cc
	misc/elsasser.cc
	misc/field_k_fns.cc
	misc/force_alpha_beta.cc
	misc/force_model.cc
	misc/getmagV0.cc
	misc/hc.cc
	misc/ic_amp_phase.cc
	misc/ic_model.cc
	misc/mult_exp_visc.cc
	misc/nusselt.cc
	misc/reality_cond.cc
	misc/sincos_horizontal.cc
	misc/zero_mode_RB_slip.cc
	EnergyTrdir/compute_cylin_ring_tr.cc  
	EnergyTrdir/compute_flux.cc
	EnergyTrdir/compute_force_feed.cc
	EnergyTrdir/compute_helicity_flux.cc
	EnergyTrdir/compute_helicity_shell_tr.cc	
	EnergyTrdir/compute_ring_tr.cc 
	EnergyTrdir/compute_shell_tr.cc	
	EnergyTrdir/energy_tr_B0.cc
	EnergyTrdir/EnergyTr.cc 
	EnergyTrdir/ET_compute_nlin.cc	
	EnergyTrdir/ET_RSprod.cc 
	EnergyTrdir/ET_transform.cc
	EnergyTrdir/fill_ring.cc
	EnergyTrdir/fill_sphere_shell.cc
	EnergyTrdir/prod_nlinV.cc 
	EnergyTrdir/skpq.cc  
	)

# TARGET_LINK_LIBRARIES(IncFlow fftw3_mpi fftw3  mpi m)

#SET(CMAKE_BUILD_TYPE distribution)
#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")
#SET(CMAKE_CXX_FLAGS "-Wall -Werror")

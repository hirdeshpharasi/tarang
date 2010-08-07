/* Tarang-4.0
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-4.0 .
 *
 * Tarang-4.0 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-4.0 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-4.0; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */


/*! \file  universal_energy.h
 * 
 * @brief  Universal functions based on basis fns
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug No known bugs
 */


#ifndef _UNIVERAL_ENERGY_H
#define _UNIVERAL_ENERGY_H

#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif

#include "universal_inline.h"


#include "../../four/four_energy.h"
#include "../../scft/scft_energy.h"


using namespace blitz;


//*********************************************************************************************

/// Total fluctuating energy of 3D complex array A
DP Get_total_energy(string basis_type, string alias_switch, int N[], Array<complx,3> A);

/// Total fluctuating energy of 3D complex arrays A and B (A.B type)
DP Get_total_energy
(
	string basis_type, 
	string alias_switch, 
	int N[],  
	Array<complx,3> A, 
	Array<complx,3> B
);;

/// Computes \f$ \sum  K^n E(k) \f$ for 3D complex array A;   k=0 excluded.
DP Get_total_Sn
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> A, 
	DP n, 
	DP kfactor[]
);

/// Compute shell spectrum Sk for 3D complex array A	
void Compute_shell_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> A, 
	DP n, 
	Array<DP,1> Sk, 
	DP kfactor[]
);

/// Compute shell spectrum Sk for 3D complex arrays A & B (A.B type)
void Compute_shell_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP n, 
	Array<DP,1> Sk, 
	DP kfactor[]
);			
			
								
/// Computes total of helicity1, helicity2, dissipation_helicity1, dissipation_helicity2.
void Compute_total_helicity
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &total_dissipation_H1, DP &total_dissipation_H2,
	DP kfactor[]
);


/// Computes shell spectrum \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$.
void Compute_shell_spectrum_helicity
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3,
	DP kfactor[]
);								

//*********************************************************************************************
// ENTROPY

/// Returns entropy of a 3D vector field A.
DP Get_entropy
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
);

/// Returns entropy of a 3D scalar field A.
DP Get_entropy_scalar(string basis_type, string alias_switch, int N[], Array<complx,3> A);


//*********************************************************************************************

/// Computes ring spectrum  of a 3D vector field A.
void Compute_ring_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);

/// Computes ring spectrum  of  3D vector fields A & B  (A.B type).									
void Compute_ring_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);

/// Computes ring spectrum  of  3D scalar field F.									
void Compute_ring_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
);		


/// real(F.G^*)
void Compute_ring_spectrum
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G,
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
);

/// Computes shell helicity spectrum \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$.
void Compute_ring_spectrum_helicity
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> H1k, 
	DP kfactor[]
);

//*********************************************************************************************

/// Computes cylinder_ring spectrum  of a 3D vector field A.
void Compute_cylinder_ring_spectrum
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);


/// Computes  cylinder_ring spectrum  of 3D vector fields A.B.
void Compute_cylinder_ring_spectrum
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);


/// Compute  cylinder_ring spectrum of 3D scalar field F.
void Compute_cylinder_ring_spectrum
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> cylinder_kpll_array,  
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
);



/// Compute  cylinder_ring spectrum of 3D scalar field A.  real(F.G^*)
void Compute_cylinder_ring_spectrum
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G, 
	Array<DP, 1> cylinder_kpll_array,  
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
);

/// Computes cylinder_ring helicity spectrum  of a 3D vector field Ai
void Compute_cylinder_ring_spectrum_helicity
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> H1k, 
	DP kfactor[]
);


//*********************************************************************************************


/// Compute spectrum imag(\f$ (\vec{B0} \cdot \vec{K'}) \vec{A}.\vec{B}^{*} \f$) 
///				for all 3D shells.
void Compute_imag_shell_spectrum_B0
(
	string basis_type,
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> Sk, 
	DP kfactor[]
);


/// Compute spectrum imag(\f$ (\vec{B0} \cdot \vec{K'}) \vec{A}.\vec{B}^{*} \f$) 
///				for all 3D rings.
void Compute_imag_ring_spectrum_B0
(
	string basis_type,
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> Sk, DP kfactor[]
);



/// Compute spectrum imag(\f$ (\vec{B0} \cdot \vec{K'}) \vec{A}.\vec{B}^{*} \f$) 
///				for all 3D cylindrical rings.
void Compute_imag_cylinder_ring_spectrum_B0
(
	string basis_type,
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> Sk, DP kfactor[]
);

#endif

//******************************** End of universal_energy.h  *********************************






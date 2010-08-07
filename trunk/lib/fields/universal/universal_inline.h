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

/*! \file  universal_inline.h
 * 
 * @brief  Universal inline functions based on basis fns
 *
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the 
 *			grid index \f$ \vec{i} \f$.
 * Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber 
 *			\f$ \vec{k} \f$ using \f$ K_i = k_i * f_i \f$ where \f$ f_i \f$ is the kfactor[i]. 
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date  Sept 2008
 * @bug   No know bugs
 */ 

 
 
#ifndef _UNIVERAL_INLINE_H
#define _UNIVERAL_INLINE_H

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

#include "../../four/four_inline.h"
#include "../../scft/scft_inline.h"

using namespace blitz;

//*********************************************************************************************

/// Get grid waveno kx given first array index l1 (local). 
inline int Get_kx(string basis_type, int l1, int N[])  
{
	if (basis_type == "FOUR")  
		return Get_kx_FOUR(l1, N);
		
	else if (basis_type == "SCFT")  
		return Get_kx_SCFT(l1, N);
		
	else 
		return 0;			// for -Wall	
} 

/// Get local array index l1 given grid waveno kx. 
inline int Get_lx(string basis_type, int kx, int N[])  
{
	if (basis_type == "FOUR")  
		return Get_lx_FOUR(kx, N);
		
	else if (basis_type == "SCFT")  
		return Get_lx_SCFT(kx, N);
		
	else 
		return 0;			// for -Wall	
}


inline int Get_ix(string basis_type, int kx, int N[])	
{
	if (basis_type == "FOUR")  
		return Get_ix_FOUR(kx, N);
	
	else if (basis_type == "SCFT")  
		return kx;
	
	else 
		return 0;			// for -Wall	
}
	
/// 3D: Get grid waveno ky given second array index l2. 	
inline int Get_ky3D(string basis_type, int l2, int N[])
{ 
	if (basis_type == "FOUR")  
		return Get_ky3D_FOUR(l2, N);
		
	else if (basis_type == "SCFT")  
		return Get_ky3D_SCFT(l2, N);
	
	else 
		return 0;			// for -Wall		
}


/// 3D: Get array index l2 given grid waveno ky.
inline int Get_ly3D(string basis_type, int ky, int N[]) 
{
	if (basis_type == "FOUR")  
		return Get_ly3D_FOUR(ky, N);
		
	else if (basis_type == "SCFT")  
		return Get_ly3D_SCFT(ky, N);
		
	else 
		return 0;			// for -Wall	
}

/// 3D: If WAVENOACTUAL \f$ |\vec{K}| \f$; else if WAVENOGRID \f$ |\vec{k}| \f$.
inline DP Kmagnitude(string basis_type, int l1, int l2, int l3, int N[], DP kfactor[])
{ 
	if (basis_type == "FOUR")  
		return Kmagnitude_FOUR(l1, l2, l3, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return Kmagnitude_SCFT(l1, l2, l3, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}

//*********************************************************************************************


/// Radius of the smallest sphere that contains the wavenumber K box Ni's. 
/// (Choices: WAVENOACTUAL or WAVENOGRID)
inline int Min_radius_outside(string basis_type, string alias_switch, int N[], DP kfactor[]) 
{
	if (basis_type == "FOUR")  
		return Min_radius_outside_FOUR(alias_switch, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return Min_radius_outside_SCFT(alias_switch, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}


/// Radius of the largest sphere that fits inside the wavenumber K box Ni's. 
/// (Choices: WAVENOACTUAL or WAVENOGRID)
inline int Max_radius_inside(string basis_type, string alias_switch, int N[], DP kfactor[]) 
{
	if (basis_type == "FOUR")  
		return Max_radius_inside_FOUR(alias_switch, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return Max_radius_inside_SCFT(alias_switch, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}


//*********************************************************************************************

/// Returns the approximate number of modes in a wavenumber K shell of radius "radius". 
///		(Choices: WAVENOACTUAL or WAVENOGRID).
inline DP Approx_number_modes_in_shell(string basis_type, int N[], int radius, DP kfactor[])
{
	if (basis_type == "FOUR")  
		return  Approx_number_modes_in_shell_FOUR(N, radius, kfactor);
		
	else if (basis_type == "SCFT")		
		return  Approx_number_modes_in_shell_SCFT(N, radius, kfactor);
		
	else 
		return 0;			// for -Wall	
}

//*********************************************************************************************

/// Modal energy in 3D
inline DP Modal_energy(string basis_type, Array<complx,3> A, int i1, int i2, int i3)
{
	if (basis_type == "FOUR")  
		return  Modal_energy_FOUR(A, i1, i2, i3);
		
	else if (basis_type == "SCFT")		
		return  Modal_energy_SCFT(A, i1, i2, i3);
		
	else 
		return 0;			// for -Wall	
}


//*********************************************************************************************
//
// FOR a given anisotropy direction.
//


/// 3D: \f$ K_{||} \f$ for anisotropy calculations.
inline DP AnisKpll(string basis_type, int i1, int i2, int i3, int N[], DP kfactor[])
{ 
	if (basis_type == "FOUR")  
		return AnisKpll_FOUR(i1, i2, i3, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return AnisKpll_SCFT(i1, i2, i3, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}


/// 3D:  \f$ K_\perp \f$ for anisotropy calculations.
inline DP AnisKperp(string basis_type, int i1, int i2, int i3, int N[], DP kfactor[])
{ 
	if (basis_type == "FOUR")  
		return AnisKperp_FOUR(i1, i2, i3, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return AnisKperp_SCFT(i1, i2, i3, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}


/// 3D:  \f$ K_{h1} \f$ for anisotropy calculations.
inline DP AnisKh1(string basis_type, int i1, int i2, int i3, int N[], DP kfactor[])
{ 
	if (basis_type == "FOUR")  
		return AnisKh1_FOUR(i1, i2, i3, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return AnisKh1_SCFT(i1, i2, i3, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}


/// 3D:  \f$ K_{h2} \f$ for anisotropy calculations.
inline DP AnisKh2(string basis_type, int i1, int i2, int i3, int N[], DP kfactor[])
{ 
	if (basis_type == "FOUR")  
		return AnisKh2_FOUR(i1, i2, i3, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return AnisKh2_SCFT(i1, i2, i3, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}

/// Cylindrical: Anis_min_Kpll
inline DP Anis_min_Kpll(string basis_type, string alias_switch, int N[], DP kfactor[])
{
	if (basis_type == "FOUR")  
		return Anis_min_Kpll_FOUR(alias_switch, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return Anis_min_Kpll_SCFT(alias_switch, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}

/// Cylindrical: Anis_max_Kpll
inline DP Anis_max_Kpll(string basis_type, string alias_switch, int N[], DP kfactor[])
{
	if (basis_type == "FOUR")  
		return Anis_max_Kpll_FOUR(alias_switch, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return Anis_max_Kpll_SCFT(alias_switch, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}

/// 3D Cylindrical: Anis_max_Krho_radius_inside the wavenumber box.
inline int Anis_max_Krho_radius_inside
(
	string basis_type, 
	string alias_switch, 
	int N[], 
	DP kfactor[]
) 			
{
	if (basis_type == "FOUR")  
		return Anis_max_Krho_radius_inside_FOUR(alias_switch, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return Anis_max_Krho_radius_inside_SCFT(alias_switch, N, kfactor);	
		
	else 
		return 0;			// for -Wall	
}


// Max polar angle
inline DP  Get_max_polar_angle(string basis_type)
{
	if (basis_type == "FOUR")  
		return Get_max_polar_angle_FOUR();
	
	else if (basis_type == "SCFT")  
		return Get_max_polar_angle_SCFT();
	
	else 
		return 0;			// for -Wall	
}



/// 3D: Returns the angle K vector makes with the anisotropic for anisotropy calculations.
inline DP AnisKvect_polar_angle
(
	string basis_type, 
	int i1, int i2, int i3, 
	int N[], 
	DP kfactor[]
)
{ 
	if (basis_type == "FOUR")  
		return AnisKvect_polar_angle_FOUR(i1, i2, i3, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return AnisKvect_polar_angle_SCFT(i1, i2, i3, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}


/// 3D: Returns the angle K vector makes with the anisotropic for anisotropy calculations.
inline DP AnisKvect_azimuthal_angle_2Din3Dgrid
(
 string basis_type, 
 int i1, int i2, int i3, 
 int N[], 
 DP kfactor[]
)
{ 
	if (basis_type == "FOUR")  
		return AnisKvect_azimuthal_angle_2Din3Dgrid_FOUR(i1, i2, i3, N, kfactor);
	
	else if (basis_type == "SCFT")  
		return AnisKvect_azimuthal_angle_2Din3Dgrid_SCFT(i1, i2, i3, N, kfactor);
	
	else 
		return 0;			// for -Wall	
}


/// 3D: Returns the azimuthal angle Kperp makes with h1 axis.
inline DP AnisKvect_azimuthal_angle
(
	string basis_type, 
	int i1, int i2, int i3, 
	int N[], 
	DP kfactor[]
)
{ 
	if (basis_type == "FOUR")  
		return AnisKvect_azimuthal_angle_FOUR(i1, i2, i3, N, kfactor);
		
	else if (basis_type == "SCFT")  
		return AnisKvect_azimuthal_angle_SCFT(i1, i2, i3, N, kfactor);
		
	else 
		return 0;			// for -Wall	
}


//*********************************************************************************************


 /*! @brief Returns the x-derivative of Vx at grid wavenumber (kx,.,.) 
  *					or \f$ \mathcal{F}(d(Vx)/dx)(kx,.,.)\f$.
  *
  *  @param  basis_type		FOUR or SCFT
  *	 @param  kfactor		factor to convert grid wavenumber to actual wavenumber
  *  @param	 kx				grid wavenumber along x
  *  @param  Vx				(*V1)(kx,.,.)
  *
  *  @return i Kx Vx in FOUR basis
  *  @return Kx Vx in SCFT basis (because Vx is Sin-transformed)
  */
inline complx DxVx(string basis_type, DP kfactor[], int kx,  complx Vx)
{
	if (basis_type == "FOUR")  
		return complex<DP>(0, kx*kfactor[1]) * Vx;
		
	else if (basis_type == "SCFT")				// Vx has odd parity
		return complex<DP>(kx*kfactor[1], 0) * Vx;
		
	else
		return complx(0.0,0.0);	
		
		
}


 /*! @brief Returns \f$ \mathcal{F}(d(Vx)/dx+d(Vy)/dy)(kx,ky,.)\f$.
  *
  *  @param  basis_type		FOUR or SCFT
  *	 @param  kfactor		factor to convert grid wavenumber to actual wavenumber
  *  @param	 kx				grid wavenumber along x
  *  @param	 ky				grid wavenumber along y
  *  @param  Vx				(*V1)(kx,ky,.)
  *  @param  Vy				(*V2)(kx,ky,.)
  *
  *  @return i (Kx Vx + Ky Vy) in FOUR basis
  *  @return (Kx Vx + i Ky Vy) in SCFT basis (because Vx is Sin-transformed)
  */
inline complx DxVx_plus_DyVy
(
	string basis_type, 
	DP kfactor[], int kx,  
	int ky, 
	complx Vx, complx Vy
)
{
	if (basis_type == "FOUR")  
		return (complex<DP>(0, kx*kfactor[1]) * Vx +  complex<DP>(0, ky*kfactor[2]) * Vy);
		
	else if (basis_type == "SCFT")				// Vx has odd parity
		return (complex<DP>(kx*kfactor[1], 0) * Vx + complex<DP>(0, ky*kfactor[2]) * Vy);
		
	else
		return complx(0.0,0.0);	  // for -Wall 	
}
	
#endif


//****************************   End of universal_inline.h  ***********************************






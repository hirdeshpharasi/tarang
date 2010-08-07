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


/*! \file scft_basic.h 
 * 
 * @brief  Contains several basis functions connected to fourier basis:  computing 
 *			modal helicity and vorticity, and A(k)/k^2 etc.
 *
 * Helicity1 = \f$ H1(K) = \vec{K} . [\vec{Vr} \times \vec{Vi}] \f$. <BR>
 * Helicity2 = \f$ H2(K) = \vec{K} . [(\vec{Vr} \times \vec{Vi})] /K^2 \f$. <BR>
 *
 * Vorticity = \f$ \Omega(K) = i \vec{K} \times \vec{V}(\vec{k}) \f$.
 *
 * \f$ \vec{V}_{Fourier} = ( -i V_x, V_y, V_z) \f$
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date 30 August 2008
 * @bug  No known bug
 */

#ifndef _SCFT_BASIC_H
#define _SCFT_BASIC_H

#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>

#include "scft_inline.h"

using namespace blitz;

//*********************************************************************************************	

/// Number of modes included inside shell(n) = \f$ [R^{sh}(n-1),  R^{sh}(n) ) \f$.
/// We also count the complex conj modes. 
int Get_local_number_modes_in_shell_SCFT(int N[], DP inner_radius, DP outer_radius, DP kfactor[]);

int Get_number_modes_in_shell_SCFT(int N[], DP inner_radius, DP outer_radius, DP kfactor[]);


/// 3D: Wavenumber components corresponding to grid index (i1, i2, i3).
void Wavenumber_SCFT(int l1, int l2, int l3, int N[], DP kfactor[], TinyVector<DP,3> &kk);


// Complex K; The imaginary part is zero.  Written to use cross function of blitz.
// Omega = cross(V,K).
void Wavenumber_SCFT(int l1, int l2, int l3, int N[], DP kfactor[], TinyVector<complx,3> &kk);

//*********************************************************************************************	

/** @brief Returns helicity of mode with grid index (i1, i2, i3)
 *
 * Given SCFT vector \f$ (V_1 , V_2, V_3) \f$, the Fourier space vector is 
 *			\f$ (V_1/i , V_2, V_3) \f$.
 * Note that the first component of the velocity vector SinFOUR transformed.  
 * We are using RB free-slip basis functions as standard for this calculation.
 *
 * @param  N[]  The size of the array A.
 * @param  l1, l2, l3 local index
 * @param  Ax  x-component of the vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  helicity = \f$ \vec{K} . [\vec{Vr} \times \vec{Vi}] \f$.
 */			
DP Get_Modal_helicity_SCFT
(
	int l1, int l2, int l3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[]
);


/** @brief Computes vorticity of mode with grid index (l1, l2, l3)
 *
 * Given SCFT vector \f$ (V_1 , V_2, V_3) \f$, the Fourier space vector is 
 *			\f$ (V_1/i , V_2, V_3) \f$.
 * Note that the first component of the velocity vector SinFOUR transformed.  
 * We are using RB free-slip basis functions as standard for this calculation. <BR>
 *
 * Vorticity is only z direction.  x and y components are zero.
 *
 * @param  N[]  The size of the array A.
 * @param  l1, l2, l3 local index
 * @param  Ax  x-component of the vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  vorticity(K) = \f$ i \vec{K} \times \vec{V}(\vec{k}) \f$.
 */											
void Compute_Modal_vorticity_SCFT
(
	int l1, int l2, int l3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[], 
	TinyVector<complx,3> &vorticity
);
								
void Compute_Modal_vorticity_y_component_SCFT
(
 int l1, int l2, int l3, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
 DP kfactor[], 
 complx &vort_y
 );

//*********************************************************************************************	

								
/// 3D: Replaces A(k) by \f$ A(k)*K^{2} \f$; A(0)=0.								
void Array_mult_ksqr_SCFT(int N[], Array<complx,3> A, DP kfactor[]);

/// 3D: Replaces A(k) by \f$ A(k)/K^{2} \f$; A(0)=0.
void Array_divide_ksqr_SCFT(int N[], Array<complx,3> A, DP kfactor[]);

/// 3D: Replaces A(k) by \f$ A(k)*\exp(factor*K^2) \f$.
void Array_exp_ksqr_SCFT(int N[], Array<complx,3> A, DP factor, DP kfactor[]);

/// 3D: Replaces A(k) by \f$ A(k)*\exp(factor*K^2 + hyperfactor*K^{4}) \f$.
void Array_exp_ksqr_SCFT(int N[], Array<complx,3> A, DP factor, DP hyper_factor, DP kfactor[]);


/// 3D: Replaces A(k) by \f$ A(k)*(\vec{V0} \cdot \vec{K})^2 / K^2 \f$.
void Array_mult_V0_khat_sqr_SCFT(int N[], Array<complx,3> A, TinyVector<DP,3> V0, DP kfactor[]);


#endif


//************************************  End of scft_basic.h ***********************************







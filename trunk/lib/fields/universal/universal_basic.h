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


/*! \file  universal_basic.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug	No known bugs
 */


#ifndef _UNIVERAL_FN_H
#define _UNIVERAL_FN_H

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif

#include <mpi.h>
#include <fftw3-mpi.h>

#include "universal_inline.h"

#include "../../four/four_basic.h"
#include "../../scft/scft_basic.h"

using namespace blitz;

//*********************************************************************************************			


/// No of modes in the shell[inner_radius, outer_radius).
int Get_number_modes_in_shell
(
	string basis_type, 
	int N[], 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
);


/// Grid wavenumber at (l1, l2, l3).
void Wavenumber
(
	string basis_type, 
	int i1, int i2, int i3, 
	int N[], 
	DP kfactor[], 
	TinyVector<DP,3> &kk
);

/// Returns helicity of mode with grid index (l1, l2, l3).
DP Get_Modal_helicity
(
	string basis_type, 
	int i1, int i2, int i3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[]
);

// Computes vorticity of mode with grid index (l1, l2, l3).
void Compute_Modal_vorticity
(
	string basis_type, 
	int i1, int i2, int i3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[], 
	TinyVector<complx,3> &vorticity
);

void Compute_Modal_vorticity_y_component
(
 string basis_type, 
 int l1, int l2, int l3, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
 DP kfactor[], 
 complx &vort_y
 );
								
								
/// 3D: Replaces A(k) by \f$ A(k)*K^{2} \f$; A(0)=0
void Array_mult_ksqr(string basis_type, int N[], Array<complx,3> A, DP kfactor[]);

/// 3D: Replaces A(k) by \f$ A(k)/K^{2} \f$; A(0)=0.
void Array_divide_ksqr(string basis_type, int N[], Array<complx,3> A, DP kfactor[]);

/// 3D: Replaces A(k) by \f$ A(k)*\exp(factor*K^2) \f$
void Array_exp_ksqr(string basis_type, int N[], Array<complx,3> A, DP factor, DP kfactor[]);

/// 3D: Replaces A(k) by \f$ A(k)*\exp(factor*K^2) \f$
void Array_exp_ksqr(string basis_type, int N[], Array<complx,3> A, DP factor, DP kfactor[]);


/// 3D: Replaces A(k) by \f$ A(k)*\exp(factor*K^2 + hyperfactor*K^{4}) \f$
void Array_exp_ksqr
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	DP factor, DP hyper_factor,
	 DP kfactor[]
);

/// 3D: Replaces A(k) by \f$ A(k)*(V0.K)^2 / K^2 \f$.
void Array_mult_V0_khat_sqr
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	TinyVector<DP,3> V0, 
	DP kfactor[]
);

																																																																																																																																																																															
								
								
#endif


//******************************** End of field_basic.h  **************************************



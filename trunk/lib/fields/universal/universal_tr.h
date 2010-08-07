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


/*! \file  universal_fn.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 */


#ifndef _UNIVERAL_TR_H
#define _UNIVERAL_TR_H

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

#include "../../four/four_tr.h"
#include "../../scft/scft_tr.h"

#include "universal_inline.h"


using namespace blitz;

//*********************************************************************************************

/// Initialization of FFTW plans for 3D real-to-complex and complex-to-real Fourier transforms.
void Init_fftw_plan_array(string basis_type, int N[], Array<complx,3> A);


//*********************************************************************************************


/*! @brief Inplace Forward transform of A.
  *
  *  @param  basis_type		FOUR or SCFT
  *  @param	 A				Real 3D array (to be transformed).
  *  @param	 temp_r			Complex 3D array, a temporary array.
  *  @param  parity			0 if Atr is to be CosFOUR transformed; 
  *							1 if Atr is to be SinFOUR transformed; 
  *
  *  @return FOUR: \f$ \mathcal{F}(A) -> A	\f$ (Inplace).
  *  @return SCFT: CFT(A) -> A if parity=0;  SFT(A) -> A if parity=1. 
  */
void Forward_transform_array
(
	string basis_type, 
	int N[],  
	Array<complx,3> A, 
	Array<complx,3> temp_r, 
	int parity
);

//
//

/*! @brief Forward transform transpose order of A.
  *
  *  @param  basis_type		FOUR or SCFT
  *  @param	 Atr			Real 3D array (to be transformed).
  *  @param  parity			0 if Atr is to be CosFOUR transformed; 
  *							1 if Atr is to be SinFOUR transformed; 
  *
  *  @return FOUR: Not defined at present.
  *  @return SCFT: \f$ \mathcal{F}(Atr) -> A \f$. A is CosFOUR transformed if parity=0.  
  *					A is SinFOUR transformed if parity=1. 
  */																		
void Forward_transform_array_transpose_order
(
	string basis_type, 
	int N[],
	Array<complx,3> Atr, 
	Array<complx,3> A, 
	int parity
);	


//*********************************************************************************************

/*! @brief Inplace Inverse transform of A.
  *
  *  @param  basis_type		FOUR or SCFT
  *  @param	 A				Complex 3D array (to be transformed).
  *  @param	 temp_r			Complex 3D array, a temporary array.
  *  @param  parity			0 if Atr is to be CosFOUR transformed; 
  *							1 if Atr is to be SinFOUR transformed; 
  *
  *  @return FOUR: \f$ \mathcal{F}^{-1}(A) -> A	\f$ (Inplace).
  *  @return SCFT: ICFT(A) -> A if parity=0;  ISFT(A) -> A if parity=1. 
  */
void Inverse_transform_array
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r, 
	int parity
);	

									

/*! @brief Inverse transform transpose order.
  *
  *  @param  basis_type		FOUR or SCFT
  *  @param	 Atr			Real 3D array (to be transformed).
  *  @param  parity			0 if A is CosFOUR transformed array; 
  *							1 if A is SinFOUR transformed array; 
  *
  *  @return FOUR: Not defined at present.
  *  @return SCFT: \f$ \mathcal{F}^{-1}(A) -> Atr	\f$.  Atr is CosFOUR transformed 
  *					if parity=0.  Atr is SinFOUR transformed if parity=1. 
  */																																																																																												
void Inverse_transform_array_transpose_order
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> Atr, 
	int parity
);


//*********************************************************************************************


/*! @brief 3D: \f$ B = \mathcal{F}(d(A)/dx) \f$.
  *
  *  @param  basis_type		FOUR or SCFT
  *	 @param  kfactor		factor to convert grid wavenumber to actual wavenumber
  *  @param	 A				Complex 3D array
  *
  *  @return B	= i Kx A	in FOUR basis
  *  @return B = \f$ (-Kx)*A \f$ in SCFT basis with parity 0 or even 
  *							(A is Cos-transformed array, B is Sin-transformed array).  
  *  @return B = \f$  Kx* A \f$ in SCFT basis with parity 1 or odd  
  *							(A is SIn-transformed array, B is Cos-transformed array).
  */
void Xderiv
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP kfactor[], 
	int parity
);	

//
//

/*! @brief 3D: \f$ B = \mathcal{F}(d(A)/dy) \f$.
  *
  *  @param  basis_type		FOUR or SCFT
  *	 @param  kfactor		factor to convert grid wavenumber to actual wavenumber
  *  @param	 A				Complex 3D array
  *
  *  @return B	= i Ky A	in FOUR/SCFT basis
  */
void Yderiv
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP kfactor[], 
	int parity
);	

//
//

/*! @brief 3D: \f$ B = \mathcal{F}(d(A)/dz) \f$.
  *
  *  @param  basis_type		FOUR or SCFT
  *	 @param  kfactor		factor to convert grid wavenumber to actual wavenumber
  *  @param	 A				Complex 3D array
  *
  *  @return B	= i Kz A	in FOUR/SCFT basis
  */
void Zderiv
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP kfactor[], 
	int parity
);	

#endif


//******************************** End of universal_tr.h  **************************************



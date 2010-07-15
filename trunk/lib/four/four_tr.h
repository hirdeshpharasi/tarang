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


/*! \file fourier.h 
 * 
 * @brief Contains declarations of Fourier transforms and derivative for arrays 
 *			in 3D for parallel code.
 * 
 * Array A is declared as complex array, but it is also used for storing real values. <BR>
 * 
  * In 3D: \f$ A[0:N_1-1, 0:N_2-1, 0:N_3/2] \f$ contains complex entries, but it could also 
 *			store \f$ A[0:N_1-1, 0:N_2-1, 0:N_3-1] \f$ real entries.  
 *			The last plane (\f$ k_z=N_3/2 \f$) is not used while storing real array. <BR>
 *
 * Since \f$ A(\vec{-k}) = A^*(\vec{k}) \f$, half the complex entries are not stored.  
 *			For example, in 3D we do not store \f$ k_z < 0 \f$.
 *
 * We use FFTW library for both forward and backward transform.  The FFTW functions perform 
 *	parallel operations.  For example,  forward Fourier transform takes all the pieces of the 
 *  arrays from all the processors and performs the transform.
 *  The full array size is inputted as parameter for FFTW functions.
 *
 * The definition of Forward transform <BR>
 * \f{equation}  A(\vec{k}) = \mathcal{F}(A(\vec{j})) 
 *			= \frac{1}{\Pi N_s} \sum_{\vec{j}} A(\vec{j}) \exp(-2 \pi i \sum j_s k_s /N_s)  \f}
 *
 * The definition of Inverse Fourier transform <BR>
 * \f{equation}  A(\vec{j}) =  \mathcal{F}^{-1}(A(\vec{k})) 
 *			=  \sum_{\vec{k}} A(\vec{k}) \exp(2 \pi i \sum j_s k_s /N_s)  \f}
 *
 * @author  M. K. Verma
 * @version 4.0 Parallel
 * @date August 2008
 * @bug		No known bugs
 */ 
 

#ifndef _FOUR_TR_H
#define _FOUR_TR_H
  
#include <blitz/array.h>
#include <complex>
#include <cmath>

#include <mpi.h>
#include <fftw3-mpi.h>

#include "four_basic.h"


using namespace blitz ;

//*********************************************************************************************

/** @brief Initialization of FFTW plans for 3D real-to-complex and 
 *			complex-to-real Fourier transforms.
 *  
 *  The 1D Fourier transforms are inplace.
 *
 * @param  A  The complex array that is Fourier transformed to real and vice versa.
 * @param  N[]  The size of the array A: \f$ (N_1, N_2, N_3) \f$.
 * 
 * @return  global variables r2c_plan_FOUR and c2r_plan_FOUR are assigned 
 *			fftw_plan_dft_r2c_3d and fftw_plan_dft_c2r_3d respectively. 
 */
void Init_fftw_plan_FOUR(int N[], Array<complx,3> A);

void Init_fftw_plan_2D_FOUR(int NN[], Array<complx,2> A);


//*********************************************************************************************

/** @brief Inplace FFTW Forward Fourier transform of 3D Array A[N].
 *
 * Since \f$ A(\vec{-k}) = A^*(\vec{k}) \f$, the negative kz values are not stored.
 * The \f$ k_x \f$ values are stored in array as \f$ 0,1, ..,N_1/2, -N_1/2+1, ...,-2,-1 \f$.
 * The \f$ k_y \f$ values are stored in array as \f$ 0,1, ..,N_2/2, -N_2/2+1, ...,-2,-1 \f$.
 * This is FFTW convention.
 *
 * @param  A  Partial Real array \f$ A [ 0:local_N1-1, 0:N_2-1, 0:N_3-1 ] \f$.
 * @param  N[]  The size of the array A:  \f$ (N_1, N_2, N_3) \f$.
 * 
 * @return  Partial Complex A(k). The full array has \f$ k_x=-N_1/2+1:N_1/2; 
 *			k_y=-N_2/2+1:N_2/2;  k_z=0:N_3/2 \f$, and it is divided over all the procs.
 */
void ArrayFFTW_FOUR(fftw_plan r2c_plan_FOUR, int N[], Array<complx,3> A); 


//*********************************************************************************************


/** @brief Inplace Forward Fourier transform of 3D Array A[N].
 *
 * Three steps: 
 * (1) Zeropad: Set the \f$ k_z = N_3/2 \f$ plane to zero, or \f$ A(:, :, N_3:N_3+1) = 0 \f$.
 * (2) FFTW forward Fourier transform.  Unnormalized.
 * (3) Normalize A(k) by dividing by \f$ N_1 * N_2 * N_3 \f$.
 *
 * @param  A  Real array \f$ A [0:local_N1-1, 0:N_2-1, 0:N_3-1 ] \f$.
 * @param  N[]  The size of the array A: \f$ (N_1, N_2, N_3) \f$.
 * 
 * @return   Partial complex A(k).   
 * 
 */
void ArrayFFT_FOUR(fftw_plan r2c_plan_FOUR, int N[], Array<complx,3> A); 

//*********************************************************************************************

/** @brief Inplace FFTW Forward Fourier transform of 3D Array A[N].
 *
 * Since \f$ A(\vec{-k}) = A^*(\vec{k}) \f$, the negative kz values are not stored.
 * The \f$ k_x \f$ values are stored in array as \f$ 0,1, ..,N_1/2, -N_1/2+1, ...,-2,-1 \f$.
 * The \f$ k_y \f$ values are stored in array as \f$ 0,1, ..,N_2/2, -N_2/2+1, ...,-2,-1 \f$.
 * This is FFTW convention.
 *
 * @param  A Complex array \f$ A [ 0:local_N1-1, -N_2/2+1:N_2/2, 0:N_3/2 ] \f$.
 * @param  N[]  The size of the array A:  \f$ (N_1, N_2, N_3) \f$.
 * 
 * @return  Partial real array \f$ A [ 0:local_N1-1, 0:N_2-1, 0:N_3-1 ] \f$.
 */
void ArrayIFFT_FOUR(fftw_plan c2r_plan_FOUR, int N[], Array<complx,3> A);


//*********************************************************************************************

/// Put zeros at the last xy-plane:  A(:,:,N3/2) = 0.
void Zero_pad_last_plane_FOUR(int N[],  Array<complx,3> A);

/// Divides A by N[1]*N[2]*N[3]
void Norm_FOUR(int N[], Array<complx,3> A);

/// Derivative along x: B(k) = i*kx*A(k)
void Xderiv_FOUR(int N[],  Array<complx,3> A, Array<complx,3> B, DP kfactor[]);

/// Derivative along y: B(k) = i*ky*A(k)
void Yderiv_FOUR(int N[], Array<complx,3> A, Array<complx,3> B, DP kfactor[]);

/// Derivative along z: B(k) = i*kz*A(k)
void Zderiv_FOUR(int N[], Array<complx,3> A, Array<complx,3> B, DP kfactor[]);

#endif

//******************************** end of four_tr.h *******************************************






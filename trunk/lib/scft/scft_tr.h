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

/*! \file sincosfour.h
 * 
 * @brief Contains declarations of SinCosFourier transforms and derivative for arrays in  3D 
 *					( MPI Implementation).
 * 
 * In 3D: \f$ A[0:N_1-1, 0:N_2-1, 0:N_3/2] \f$ contains complex entries, but it could also 
 *			store \f$ A[0:N_1-1, 0:N_2-1, 0:N_3-1] \f$ real entries.  The last plane 
 *			(\f$ k_z=N_3/2 \f$) is not used while storing real array. <BR>
 *
 * We use FFTW library for both forward and backward transform.
 *
 *  1D Forward Sin transform (RODFT10 or SFT-II)   <BR>
 * \f{equation} A(k) =  \mathcal{F}(A(j)) =  \frac{1}{2 N_1} \sum_{j=0}^{N_1-1} A(j) 
 *											2 \sin [ \pi (j+1/2) k /N_1 ]  \f}
 *
 * 1D Inverse Sin transform (RODFT01 or SFT-III) <BR>
 * \f{equation} A(j) =  \mathcal{F}(A(k)) =   \sum_{m=1}^{N_1-1} A(k)  
 *												2 \sin [ \pi (j+1/2) k /N_1 ] \f}
 *
 *
 * 1D Forward Cos transform (REDFT10 or DCT-II) <BR>
 * \f{equation}  A(k) =  \mathcal{F}(A(j)) =  \frac{1}{2 N_1} \sum_{j=0}^{N_1-1} A(j) 
 *												2 \cos [ \pi (j+1/2) k /N_1 ]  \f}
 *
 * 1D Inverse Cos transform (REDFT01 or DCT-III) <BR>
 * \f{equation}   A(j) =  \mathcal{F}(A(k)) =  A(0) +  \sum_{j=0}^{N_1-1} A(k) 
 *														2 \cos [ \pi (j+1/2) k /N_1 ]   \f}
 *
 * Forward SinFourier transform  (s sum along perp dirns) <BR>
 * \f{equation}  A(\vec{k}) = \frac{1}{2 N_1 \Pi  N_s} \mathcal{F}(A(\vec{j})) =  
 *			\sum_{\vec{k}} A(\vec{j}) 2\sin [ \pi (j_{||}+1/2) k_{||} /N_1 ]  
 *				\exp(2 \pi i \sum j_s k_s /N_s) \nonumber  \f} 
 *
 * The definition of Inverse SinFourier transform (s sum along perp dirns) <BR>
 * \f{equation}  A(\vec{j}) = \mathcal{F}^{-1}(A(\vec{k})) =  \sum_{\vec{k}} A(\vec{k}) 
 *			2\sin [ \pi (j_{||}+1/2) k_{||} /N_1 ] 
 *			\exp(2 \pi i \sum j_s k_s /N_s) \nonumber  \f}  
 *
 * Forward CosFourier transform  (s sum along perp dirns) <BR>
 * \f{equation}  A(\vec{k}) = \frac{1}{2 N_1 \Pi N_s} \mathcal{F}(A(\vec{j})) 
 *					=  \sum_{\vec{k}} A(\vec{j}) 2\cos [ \pi (j_{||}+1/2) k_{||} /N_1 ]  
 *						\exp(2 \pi i \sum j_s k_s /N_s) \nonumber  \f} 
 *
 * The definition of Inverse CosFourier transform (s sum along perp dirns) <BR>
 * \f{equation}  A(\vec{j}) = \mathcal{F}^{-1}(A(\vec{k})) =  \sum_{\vec{k}_{\perp}} 
 *								A(0,\vec{k}_{\perp}) \exp(2 \pi i \sum j_s k_s /N_s)  
 *							 + \sum_{\vec{k}} A(\vec{k}) 2\cos [ \pi (j_{||}+1/2) k_{||} /N_1 ] 
 *								\exp(2 \pi i \sum j_s k_s /N_s) \nonumber  \f} 
 * 
 *	1D transforms are real to real. The range of j and k = [0..N1-1]. <BR>
 *  2D and 3D forward transforms are real to complex, and the inverse transforms are 
 *			from complex to real. <BR>
 *
 *  FFTW's and our definitions differ a bit for the Sin transform. <BR>
 *  For forward Sin transform, FFTW saves k = 1..n in the array; for forward Cos transform, 
 *			FFTW saves k = 0..n-1.  
 *  To keep consistency, we save k = 0..n-1 for both Sin/Cos transform.  
 *			For Sin transform k=0 mode is zero.
 *  Due to the above reason, for FFTW-Sin to Our Sin transform --->  Shift right.  <BR>
 *  for Our Sin transform to FFT-Sin transform --> Shift left.
 *
 * The reality condition implies that \f$ A(k_{||}, \vec{-k}_{\perp}) 
 *			= (A(k_{||})^*, \vec{k}_{\perp}) \f$. Hence only half modes are stored. In 2D, 
 *				ky >= 0, while in 3D, kz >= 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date	August 2008
 * @bug		No known bugs
 */ 

#ifndef _SCFT_TR_H
#define _SCFT_TR_H
  
#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>

#include "scft_basic.h"


using namespace blitz ;


//*********************************************************************************************

/** @brief Initialization of FFTW plans for 3D Sin/Cos FOUR forward and inverse transforms.
 *  
 *  The 2D Fourier transforms are inplace.  Sin/Cos along x direction, and FOUR along yz dirns.
 *
 * @param  A  The complex array that is Sin/CosFourier transformed to real and vice versa.
 * @param  N[]  The size of the array A.
 * 
 * @return  global variables sintr_plan_SCFT = fftw_plan_dft_r2c_1d (FFTW_RODFT10) 
 *												-- for x direction. 
 * @return  global variables costr_plan_SCFT = fftw_plan_dft_r2c_1d (FFTW_REDFT10) 
 *												-- for x direction. 
 * @return  global variables isintr_plan_SCFT = fftw_plan_dft_r2c_1d (FFTW_RODFT01) 
 *												-- for x direction.
 * @return  global variables isintr_plan_SCFT = fftw_plan_dft_r2c_1d (FFTW_REDFT01) 
 *												-- for x direction.
 * @return  global variables r2c_plan_FOUR = fftw_plan_dft_r2c_2d -- for yz direction. 
 * @return  global variables c2r_plan_FOUR = fftw_plan_dft_c2r_2d -- for yz direction. 
 */
void Init_fftw_plan_SCFT(int NN[], Array<complx,3> A);			


//*********************************************************************************************

/// Put zeros at the last xy-plane:  A(:,:,N3/2) = 0.
void Zero_pad_lastplane_SCFT(int N[], Array<complx,3> A);			


/// Divides A by 2*N[1]*N[2]*N[3].
void Norm_SCFT(int N[], Array<complx,3> A);						

/// 1D FFTW-Sin transform Row().
void Sintr_row_SCFT(fftw_plan sintr_plan_SCFT, int N[], Array<DP,1> Row);

/// 1D FFTW-Cos transform Row().
void Costr_row_SCFT(fftw_plan costr_plan_SCFT, int N[], Array<DP,1> Row);

/// 1D FFTW-InverseSin transform Row().
void ISintr_row_SCFT(fftw_plan isintr_plan_SCFT, int N[], Array<DP,1> Row);

/// 1D FFTW-InverseCos transform Row().
void ICostr_row_SCFT(fftw_plan icostr_plan_SCFT, int N[], Array<DP,1> Row);

/// 2D FFTW Forward Fourier transform Plane().
void FT_Plane_SCFT(fftw_plan r2c_plan_SCFT, int N[], Array<complx,2> Plane);

/// 2D FFTW Inverse Fourier transform Plane().	
void IFT_Plane_SCFT(fftw_plan c2r_plan_SCFT, int N[], Array<complx,2> Plane);

void FT_column_1d_SCFT(fftw_plan r2c_1d_plan_SCFT, int N[], Array<complx,1> temp_column_2d);

void IFT_column_1d_SCFT(fftw_plan c2r_1d_plan_SCFT, int N[], Array<complx,1> temp_column_2d);

/// Shift-right array A.  Before the shift kx=1..n.
/// After the shift: kx=0..(n-1).  A(0,:,:)=0. kx=n mode thrown out.						 	 
void ArrayShiftRight_SCFT(int N[], Array<complx,3> A);

/// Shift-left array A.  Before the shift kx=0..(n-1).
/// After the shift: kx=1..(n).  A(n,:,:)=0. 
void ArrayShiftLeft_SCFT(int N[], Array<complx,3> A);

//*********************************************************************************************

/** @brief SinFour transform: SFT(Atr) = A.
 *
 * Steps: <BR>
 * (1) FFTW forward Sin transform Atr along y directions for all columns.  Unnormalized. ky=1..n. <BR>
 * (2) Zero_pad last plane of Atr (k = N[3]/2). <BR>
 * (3) Inverse transpose Atr, and place the result in A. Transpose along x<->y. <BR>
 * (4) Fourier transform of all yz planes. <BR>
 * (5) Normalize A(k) by dividing by \f$ 2 N_1 N_2 N_3 \f$.  <BR>
 * (6) Shift right. After the shift kx=0..(n-1).
 *
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$. 
 * @param  Atr(j)  Real array \f$ Atr[0:localN2-1, 0:N_1-1, 0:N_3-1] \f$. The complete array  \f$ Atr [0:localN2-1, 0:N_1-1, 0:N_3-1] \f$.
 * 
 * @return A(k) with \f$ k_x=0:local_N1-1, k_y=0:N_2-1, k_z=0:N-3/2 \f$.  
 */
void ArraySFT_SCFT_transpose_order
(
	fftw_plan sintr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT, 
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> Atr, 
	Array<complx,3> A
);

//*********************************************************************************************
											
/** @brief CosFour transform  CFT(Atr) = A.
 *
 * Steps: <BR>
 * (1) FFTW forward Cos transform along y directions for all columns.  Unnormalized. kx=0..n-1. <BR>
 * (2) Zero_pad last plane of Atr (k = N[3]/2). <BR>
 * (3) Inverse transpose Atr, and place the result in A. Transpose along x<->y. <BR>
 * (4) Fourier transform of all yz planes. <BR>
 * (4) Normalize A(k) by dividing by \f$ 2 N_1 N_2 N_3\f$.  
 *
 * @param  Atr(j)  Real array \f$ A [0:N_1-1, 0:N_2-1, 0:N_3-1] \f$.
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$.
 * 
 * @return A(k) with \f$ k_x=0:N_1-1,k_y=0:N_2-1, k_z=0:N-3/2 \f$.  
 */											
void ArrayCFT_SCFT_transpose_order
(
	fftw_plan costr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT,
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> Atr, 
	Array<complx,3> A
);

//*********************************************************************************************

/** @brief InverseSinFour transform ISFT(A) -> Atr.
 *
 * Steps: <BR>
 * (1) Shift left. After the shift kx=1..(n). <BR>
 * (2) Inverse Fourier transform of all yz planes. <BR>
 * (3) Transpose(A) -> Atr.
 * (4) FFTW Inverse Sin transform along y directions for all columns.  i1=0..(n-1). <BR>
 *
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$. 
 * @param  A(k)  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$.
 * 
 * @return \f$ Atr(0:localN2-1, 0:N_1-1, 0:N_3-1) \f$ real array.  
 */
void ArrayISFT_SCFT_transpose_order
(
	fftw_plan isintr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT, 
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> Atr
);

//*********************************************************************************************	
																					
/** @brief InverseCosFour transform ICFT(A) -> Atr.
 *
 * Steps: <BR>
 * (1) Inverse Fourier transform of all yz planes. <BR>
 * (2) Transpose(A) -> Atr.
 * (3) FFTW Inverse Cos transform along y directions for all columns.  i1=0..(n-1). <BR>
 *
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$. 
 * @param  A(k)  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3/2] \f$.
 * 
 * @return \f$ Atr(0:localN2-1, 0:N_1-1, 0:N_3-1) \f$ real array.  
 */									
void ArrayICFT_SCFT_transpose_order
(
	fftw_plan icostr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT, 
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> Atr
); 

//*********************************************************************************************
											
//
// IN PLACE
//

/** @brief Inplace Forward SinFOUR transform of 3D real Array A[].
 *
 * Steps: <BR>
 * (1) Transpose(A) -> temp_r <BR>
 * (2) ArraySFT_SCFT_transpose_order(sintr_plan_SCFT, r2c_plan_SCFT, N, temp_r, A). 
 *		That is, SFT(temp_r) -> A.
 *
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$.
 * @param  A(j)  Real array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$.
 * @param  temp_r	temporary array with dim [localN2, N1, N3/2+1].
 *
 * @return A(k): SFT of the original array. \f$ k_x=0:localN1-1, k_y=0:N_2-1, k_z=0:N-3/2 \f$.  
 */
void ArraySFT_SCFT
(
	fftw_plan sintr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT, 
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
);

//*********************************************************************************************
										
/** @brief Inplace Forward CosFOUR transform of 3D real Array A[].
 *
 * Steps: <BR>
 * (1) Transpose(A) -> temp_r <BR>
 * (2) ArrayCFT_SCFT_transpose_order(costr_plan_SCFT, r2c_plan_SCFT, N, temp_r, A). 
 *		That is, CFT(temp_r) -> A.
 *
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$.
 * @param  A(j)  Real array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$.
 * @param  temp_r	temporary array with dim [localN2, N1, N3/2+1].
 *
 * @return A(k): CFT of the original array. \f$ k_x=0:localN1-1, k_y=0:N_2-1, k_z=0:N-3/2 \f$.  
 */									 
void ArrayCFT_SCFT
(
	fftw_plan costr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT, 
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
);		

//*********************************************************************************************

									
/** @brief Inplace Inverse SinFOUR transform of 3D complex Array A[].
 *
 * Steps: <BR>
 * (1) ArrayISFT_SCFT_transpose_order(sintr_plan_SCFT, r2c_plan_SCFT, N, A, temp_r). 
 *		That is, ISFT(A) -> temp_r.
 * (2) Transpose(temp_r) -> A.
 *
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$.
 * @param  A(k)  complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$.
 * @param  temp_r	temporary array with dim [localN2, N1, N3/2+1].
 *
 * @return A(j): A real array that is ISFT of the original array. 
 *				\f$ k_x=0:localN1-1, k_y=0:N_2-1, k_z=0:N-3/2 \f$.  
 */																																																																
void ArrayISFT_SCFT
(
	fftw_plan isintr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT, 
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
);


//*********************************************************************************************
										
/** @brief Inplace Inverse CosFOUR transform of 3D complex Array A[].
 *
 * Steps: <BR>
 * (1) ArrayICFT_SCFT_transpose_order(costr_plan_SCFT, r2c_plan_SCFT, N, A, temp_r). 
 *		That is, ISFT(A) -> temp_r.
 * (2) Transpose(temp_r) -> A.
 *
 * @param  N[]  The size of the array A: \f$ (N_1, N_2) \f$.
 * @param  A(k)  complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N_3-1] \f$.
 * @param  temp_r	temporary array with dim [localN2, N1, N3/2+1].
 *
 * @return A(j): A real array that is ISFT of the original array. 
 *					\f$ k_x=0:localN1-1, k_y=0:N_2-1, k_z=0:N-3/2 \f$.  
 */										
void ArrayICFT_SCFT
(
	fftw_plan icostr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT, 
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
);
																				


//*********************************************************************************************


/** @brief  Derivative along x of 3D SinFOUR-transformed array.
 *
 * @param  A  Complex array \f$ A [0:loalN1-1, 0:N_2-1, 0:N-3/2] \f$.
 * @param  N[]  The size of the array A: \f$ (N_1, N_2, N_3) \f$.
 * 
 * @return B = \f$ f_1 * i1 * A \f$  (B is Cos-transformed array).  
 */										
void Xderiv_Sin_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[]); 

//*********************************************************************************************

/** @brief  Derivative along x of 3D CosFOUR-transformed array.
 *
 * @param  A  Complex array \f$ A [0:localN1-1, 0:N_2-1, 0:N-3/2] \f$.
 * @param  N[]  The size of the array A: \f$ (N_1, N_2, N_3) \f$.
 * 
 * @return B = \f$ -f_1 * i1 * A \f$  (B isSin-transformed array).  
 */
void Xderiv_Cos_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[]); 

//*********************************************************************************************

/// 3D: Derivative along y: B(k) = i*Ky*A(k)
void Yderiv_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[]);

//*********************************************************************************************


/// 3D: Derivative along z: B(k) = i*Kz*A(k)
void Zderiv_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[]); 


#endif
 
//******************************** end of scft_tr.h *******************************************




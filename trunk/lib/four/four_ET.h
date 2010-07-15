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



/*! \file four_ET.h 
 * 
 * @brief Contains useful functions for computing products of vector fields over shells and rings.
 * 
 *  Data array is complex A:  \f$ A[0:N_1-1, 0:N_2-1, 0:N_3/2] \f$. <BR>
 *
 *  We compute sum  \f$  \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) /2 ] \f$  for over a region.  
 *  This region could be a shell, ring or a cylinderial segment. <BR>
 *
 *  We compute sum  \f$  \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) ] \f$  for over a region.  
 * This region could be a shell, ring or a cylinderial segment. <BR>
 *
 *	shell(n) = \f$  K' \in (R^{sh}(n-1), R^{sh}(n) ] \f$.  
 *	with typical \f$ R^{sh} = 0,2,4,8,..., R_{max},\infty \f$.
 *
 *	ring(n,m) = \f$ K' = (R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ] \f$.
 *	with typical \f$ R^{sh} = 0,4,8,..., R_{max},\infty \f$ and
 *		\f$ \Theta = [0:\Theta_{max}] \f$ divided in m segments.
 *
 * @sa four_inline.h
 * @sa four_tr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * 
 * @bugs  No known bug
 */ 
 
 

#ifndef _FOUR_ET_H
#define _FOUR_ET_H

#include <blitz/array.h>
#include <complex>
#include <cmath>

#include <mpi.h>
#include <fftw3-mpi.h>

#include "four_basic.h"

using namespace blitz;
				

//*********************************************************************************************

/** @brief Compute local product Real(A.B*) for a 3D shell = [inner_radius, outer_radius).
 *  
 * For spectrum calculation (\f$  Real(A.B^*)/2 \f$) we multiply the terms 
 *	by Multiplicity_factor_FOUR(). Hence, for shell_mult (\f$  Real(A.B^*) \f$) we need to 
 *  multiply the terms by 2*Multiplicity_factor_FOUR().  
 *	Modes on the surface of inner_radius is included in the shell.
 *	Exclude origin when multiplying the first shell.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * @param  B  a complex array.
 * @param  inner_radius  The inner radius of the shell.
 * @param  outer_radius  The outer radius of the shell.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 * 
 * @return  local_prod_sum = \f$  \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) ] \f$ 
 *			for \f$ K' \in \f$ (inner_radius, outer_radius].
 */								
DP Local_shell_mult_single_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
);
	
															
/// Sum local shell-mult results and place it in the master node.									
DP Shell_mult_single_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
);								

//*********************************************************************************************

/** @brief Compute local products Real(A).Real(B) and Imag(A).Imag(B) for a 3D shell
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * @param  B  a complex array.
 * @param  inner_radius  The inner radius of the shell.
 * @param  outer_radius  The outer radius of the shell.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 * 
 * @return  local_real =  \f$  \sum[ \Re(A(\vec{K}')) * \Re(B^{*}(\vec{K}')) ] \f$   
 * @return  local_imag = \f$ \sum[ \Im(A(\vec{K}')) * \Im(B^{*}(\vec{K}')) ] \f$.
 *
 */								
void Local_shell_mult_single_real_imag_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP& local_real, DP& local_imag,
	DP kfactor[] 
);


/** @brief Compute the total products Real(A).Real(B) and Imag(A).Imag(B) for a 3D shell  
 *				by summing the local_prod_sum.
 *
 * @return  total_real =  \f$  \sum[ \Re(A(\vec{K}')) * \Re(B^{*}(\vec{K}')) /2 ] \f$   
 * @return  total_imag = \f$ \sum[ \Im(A(\vec{K}')) * \Im(B^{*}(\vec{K}')) /2 ] \f$ by 
 *							summing local_real(n) and local_imag(n).
 *
 * @sa Local_shell_mult_single_real_imag_FOUR(int N[],  Array<complx,3> A, Array<complx,3> B, 
 *			DP inner_radius, DP outer_radius, DP kfactor[], DP& local_real, DP& local_imag)
 */	
void Shell_mult_single_real_imag_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP& total_real, DP& total_imag,
	DP kfactor[]
);


//*********************************************************************************************

/** @brief Compute local product Real(A.B*) for all 3D shell.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * @param  B  a complex array.
 * @param  shell_radius_array The shell radii are stored here 
 *						\f$ R^{sh} = 0,2,4...,R_{max},\infty \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 * 
 * @return  local_result(n) = \f$  \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) /2 ] \f$ for 
 *				\f$ K' \in [R^{sh}(n-1), R^{sh}(n) ) \f$,
 */	
void Local_shell_mult_all_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result, 
	DP kfactor[]
);

/** @brief Compute product Real(A.B*) for all 3D shells by summing local_result(n).
 * 
 * @return  result(n) = \f$  \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) /2 ] \f$ 
 *				by summing local_result(n) .
 *
 * @sa Shell_mult_single_FOUR(int N[], Array<complx,2> A, Array<complx,2> B, 
 *						DP inner_radius, DP outer_radius, DP kfactor[])
 */
void Shell_mult_all_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
);


//*********************************************************************************************

/** @brief Compute local products Real(A).Real(B) and Imag(A).Imag(B) for all 3D shells by 
*					summing local_result(n).
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * @param  B  a complex array.
 * @param  shell_radius_array The shell radii are stored here 
 *					\f$ R^{sh} = 0,2,4...,R_{max},\infty \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 * 
 * @return  local_result_real(n) =  \f$  \sum[ \Re(A(\vec{K}')) * \Re(B^{*}(\vec{K}')) /2 ] \f$   
 * @return  local_result_imag(n) = \f$ \sum[ \Im(A(\vec{K}')) * \Im(B^{*}(\vec{K}')) /2 ] \f$ 
 *	for \f$ K' \in [R^{sh}(n-1), R^{sh}(n) ) \f$,
 */	
void Local_shell_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result_real, Array<DP,1> local_result_imag, 
	DP kfactor[]
);

/** @brief Compute products Real(A).Real(B) and Imag(A).Imag(B) for all 3D shell.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * @param  B  a complex array.
 * @param  shell_radius_array The shell radii are stored here 
 *				\f$ R^{sh} = 0,2,4...,R_{max},\infty \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 * 
 * @return  result_real(n) =  \f$  \sum[ \Re(A(\vec{K}')) * \Re(B^{*}(\vec{K}')) /2 ] \f$   
 * @return  result_imag(n) = \f$ \sum[ \Im(A(\vec{K}')) * \Im(B^{*}(\vec{K}')) /2 ] \f$ 
 *			by summing up local_result_real(n) and local_result_imag(n).
 *
 */	
void Shell_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result_real, Array<DP,1> result_imag, 
	DP kfactor[]
);


//*********************************************************************************************

/** @brief Compute local products \f$ 2* \Re( \vec{A}.\vec{B}^{*} )\f$ over all 3D shells.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bx  x-component of the vector field \f$ \vec{B} \f$.
 * @param  By  y-component of the vector field \f$ \vec{B} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{B} \f$.
 * @param  shell_radius_array The shell radii are stored here 
 *						\f$ R^{sh} = 0,2,4...,R_{max},\infty \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 *
 * @return  local_result(n) = \f$   \sum[ \Re(\vec{A}(\vec{K}') 
 *			\cdot \vec{B}^{*}(\vec{K}')) /2 ] \f$ for \f$ K' \in [R^{sh}(n-1), R^{sh}(n) ) \f$,
 *
 */		
void Local_shell_mult_all_FOUR
(	
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result, 
	DP kfactor[]
);



/** @brief Compute local products \f$ 2* \Re( \vec{A}.\vec{B}^{*} )\f$ for all 3D 
 *		shells by summing local_result(n).
 *
 * @return  result(n) = \f$   \sum[ \Re(\vec{A}(\vec{K}') \cdot \vec{B}^{*}(\vec{K}')) /2 ] \f$ 
 *			by summing up local_result(n).
 *
 */	
void Shell_mult_all_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,   
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
);

//*********************************************************************************************	

/** @brief Compute local products real(\f$ \vec{A} \f$).real(\f$ \vec{B} \f$) 
 *			and imag(\f$ \vec{A} \f$).imag(\f$ \vec{B} \f$) over all 3D shells.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bx  x-component of the vector field \f$ \vec{B} \f$.
 * @param  By  y-component of the vector field \f$ \vec{B} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{B} \f$.
 * @param  shell_radius_array The shell radii are stored here 
 *				\f$ R^{sh} = 0,2,4...,R_{max},\infty \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 *
 * @return local_result_real(n) = \f$ \sum[ \Re(\vec{A}(\vec{K}')) 
 *									* \Re(\vec{B}^{*}(\vec{K}')) /2 ] \f$  
 * @return local_result_imag(n) = \f$ \sum[ \Im(\vec{A}(\vec{K}')) 
 *									* \Im(\vec{B}^{*}(\vec{K}')) /2 ] \f$, 
 *									for \f$ K' \in [R^{sh}(n-1), R^{sh}(n) ) \f$,
 */	
void Local_shell_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,   
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result_real, Array<DP,1> local_result_imag, 
	DP kfactor[]
);


/** @brief Compute products \f$ 2 \Re(\vec{A}).\Re(\vec{B}) \f$ 
 *					and \f$ 2 \Im(\vec{A} ).\Im( \vec{B}) \f$
 *					for all 3D shells by summing local_result(n).
 *
 * @return result_real(n) = \f$\sum[ \Re(\vec{A}(\vec{K}')) * \Re(\vec{B}^{*}(\vec{K}')) /2 ]\f$  
 * @return result_imag(n) = \f$\sum[ \Im(\vec{A}(\vec{K}')) * \Im(\vec{B}^{*}(\vec{K}')) /2 ]\f$  
 *			by summing up local_result_real(n) and local_result_imag(n).
 *
 */
void Shell_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result_real, Array<DP,1> result_imag, 
	DP kfactor[]
);

/**********************************************************************************************

			RING MULT

***********************************************************************************************/


/*! @brief Compute local products real(A.B*) over all 3D rings.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * @param  B  a complex array.
 * @param  ring_shell_radius_array The shell radii of the rings are stored here: 
 *			\f$ R^{sh} = 0,4...,R_{max},\infty \f$.
 * @param  sector_angle_array The sector angles are stored here: 
 *			\f$ \Theta = 0,\Theta_{max}\f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 * 
 * @return  local_result(n,m) = \f$ \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) /2 ] \f$ for 
 *				\f$ K' = [R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ) \f$.
 *
 */
void Local_ring_mult_all_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
);																															


/*! @brief Compute products real(A.B*) for all 3D rings by summing local_result(n).
 *
 * @return  result(n) = \f$   \sum[ \Re(A(\vec{K}') \B^{*}(\vec{K}')) /2 ] \f$ 
 *			by summing up local_result(n).
 *
 */	
void Ring_mult_all_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
);

//*********************************************************************************************


/** @brief Compute local products Real(A).Real(B) and Real(A).Real(B) over all 3D rings.
 *
 * @return  local_result_real(n,m) = \f$  \sum[ \Re(A(\vec{K}')) *\Re(B^{*}(\vec{K}')) /2] \f$  
 * @return  local_result_imag(n,m) = \f$  \sum[ \Im(A(\vec{K}')) * \Im(B^{*}(\vec{K}'))/2]\f$  
 *				for \f$ K' = [R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ) \f$.
 *
 */	
void Local_ring_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
);

/** @brief Compute products Real(A).Real(B) and Real(A).Real(B) over all 3D rings by 
 *			summing local_result(n).
 *
 * @return  result_real(n,m) =  \f$  \sum[ \Re(A(\vec{K}')) * \Re(B^{*}(\vec{K}')) /2 ] \f$  
 * @return  result_imag(n,m) =  \f$  \sum[ \Im(A(\vec{K}')) * \Im(B^{*}(\vec{K}')) /2 ] \f$ for  
 *					by summing up local_result_real(n) and local_result_imag(n).
 *
 */	
void Ring_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
);								

//*********************************************************************************************

/** @brief Compute local products real(\f$ \vec{A}.\vec{B}^{*} \f$) over all 3D rings.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bx  x-component of the vector field \f$ \vec{B} \f$.
 * @param  By  y-component of the vector field \f$ \vec{B} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{B} \f$.
 * @param  ring_shell_radius_array The shell radii of the rings are stored here: 
 *				\f$ R^{sh} = 0,4...,R_{max},\infty \f$.
 * @param  sector_angle_array The sector angles are stored here: 
 *				\f$ \Theta = 0: \Theta_{max} \f$ in $m$ intervals.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return local_result(n,m) =\f$  \sum[ \Re(\vec{A}(\vec{K}') * \vec{B}^{*}(\vec{K}')) /2 ]\f$  
 *			for \f$ K' = [R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ) \f$.
 */		
void Local_ring_mult_all_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
);

/** @brief Compute products real(\f$ \vec{A}.\vec{B}^{*} \f$) over all 3D rings by 
 *			summing local_result(n).
 *
 * @return  result(n,m) =\f$  \sum[ \Re(\vec{A}(\vec{K}') * \vec{B}^{*}(\vec{K}')) /2 ] \f$ for 
 *			by summing up local_result_real(n).
 *
 * 
 */	
void Ring_mult_all_FOUR
(
	string alias_switch,
	int N[],  
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
);

//*********************************************************************************************


/** @brief Compute local products real(\f$ \vec{A} \f$).real(\f$ \vec{B} \f$) 
 *			and imag(\f$ \vec{A} \f$).imag(\f$ \vec{B} \f$) over all 3D rings.
 *
 * @return  result_real(n) =\f$ \sum[ \Re(\vec{A}(\vec{K}')) * \Re(\vec{B}^{*}(\vec{K}')) /2]\f$  
 * @return	result_imag(n) =\f$ \sum[ \Im(\vec{A}(\vec{K}')) * \Im(\vec{B}^{*}(\vec{K}')) /2]\f$    
 *				for \f$ K' = [R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ) \f$.
 *
 */ 
void Local_ring_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
);								


/** @brief Compute products real(\f$ \vec{A} \f$).real(\f$ \vec{B} \f$) 
 *			and imag(\f$ \vec{A} \f$).imag(\f$ \vec{B} \f$) over all 3D rings
 *			by summing local_result(n).
 *
 * @return  result_real(n) =\f$ \sum[ \Re(\vec{A}(\vec{K}')) * \Re(\vec{B}^{*}(\vec{K}')) /2]\f$  
 * @return	result_imag(n) =\f$ \sum[ \Im(\vec{A}(\vec{K}')) * \Im(\vec{B}^{*}(\vec{K}')) /2]\f$      
 *				by summing up local_result_real(n) and local_result_imag(n).
 *
 */ 
void Ring_mult_all_real_imag_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
);

/**********************************************************************************************

			Cylinderical rings
			
***********************************************************************************************/


/** @brief Compute products real(A.B*) over all 3D cylindrical rings.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * @param  B  a complex array.
 * @param  ring_shell_radius_array The shell radii of the rings are stored here: 
 *			\f$ R^{sh} = 0,4...,R_{max},\infty \f$.
 * @param  cylinder_kpll_array: [-min_kpll, ..., max_kpll].
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 * 
 * @return result(n,m) = \f$  \sum[ \Re(A(\vec{K}') * B^{*}(\vec{K}')) ] \f$ for 
 *					\f$ K' = (R^{ring}(n-1), H(m-1); R^{ring}(n), H(m) ] \f$.
 *
 */ 
void Local_cyl_ring_mult_all_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
);																								
							
	
/// Sum local ring sums and send to the master node.
void Cyl_ring_mult_all_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
);


//*********************************************************************************************

/** @brief Compute products Real(A).Real(B) and Real(A).Real(B) over all 3D rings.
 *
 * @return  result_real(n,m) =  \f$  \sum[ \Re(A(\vec{K}')) * \Re(B^{*}(\vec{K}'))  ] \f$  
 * @return  result_imag(n,m) = 	\f$  \sum[ \Im(A(\vec{K}')) * \Im(B^{*}(\vec{K}'))  ] \f$ for  
 *			\f$ K' = (R^{ring}(n-1), H(m-1); R^{ring}(n), H(m) ] \f$.
 *
 */	
void Local_cyl_ring_mult_all_real_imag_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
);


/// Sum local ring sums and send to the master node.
void Cyl_ring_mult_all_real_imag_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
);


//*********************************************************************************************


/** @brief Compute products real(\f$ \vec{A}.\vec{B}^{*} \f$) over all 3D rings.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bx  x-component of the vector field \f$ \vec{B} \f$.
 * @param  By  y-component of the vector field \f$ \vec{B} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{B} \f$.
 * @param  shell_radius_array The shell radii are stored here.
 * @param  ring_shell_radius_array The shell radii of the rings are stored here.
 * @param  sector_angle_array The sector angles are stored here:
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  result(n,m) =\f$  \sum[ \Re(\vec{A}(\vec{K}') * \vec{B}^{*}(\vec{K}')) ] \f$ for 
 *				\f$ K' = (R^{ring}(n-1), H(m-1); R^{ring}(n), H(m) ] \f$.
 */	
void Local_cyl_ring_mult_all_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
);


/// Sum local ring sums and send to the master node.
void Cyl_ring_mult_all_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
);

//*********************************************************************************************

/** @brief Compute products real(\f$ \vec{A} \f$).real(\f$ \vec{B} \f$) and 
 *						imag(\f$ \vec{A} \f$).imag(\f$ \vec{B} \f$) over all 3D cyl_rings.
 *
 * @return  result_real(n) = \f$ \sum[ \Re(\vec{A}(\vec{K}')) * 
 *									\Re(\vec{B}^{*}(\vec{K}'))  ] \f$  
 * @return	result_imag(n) = \f$  \sum[ \Im(\vec{A}(\vec{K}')) * 
 *									\Im(\vec{B}^{*}(\vec{K}')) ] \f$ for   
 *					\f$ K' = (R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ] \f$.
 *
 */ 
void Local_cyl_ring_mult_all_real_imag_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> local_result_real, Array<DP,2> local_result_imag, 
	DP kfactor[]
);


/// Sum local ring and send to the master node.
void Cyl_ring_mult_all_real_imag_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array, 
	Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
);


/**********************************************************************************************

			(B0.k) Im[vectA. conj(vectB)]

***********************************************************************************************/


/** @brief Compute \f$ (\vec{B_0} \cdot \vec{K}) \Im [\vec{A}(K) \cdot \conj{\vec{B}(K) ] \f$
 *			over a 3D shell.
 *
 * @param  N[]  The size of the array A.
 * @param  B0	The mean magnetic field.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bx  x-component of the vector field \f$ \vec{B} \f$.
 * @param  By  y-component of the vector field \f$ \vec{B} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{B} \f$.
 * @param  shell_radius_array The shell radii are stored here.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  result(n) = \f$ (\vec{B_0} \cdot \vec{K}) 
 *							\Im [\vec{A}(K) \cdot \conj{\vec{B}(K) ] \f$  for shell n.
 */ 
void Local_shell_mult_all_imagVW_B0_FOUR
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> local_result, 
	DP kfactor[]
);

/// Sum over local shells and send to the master node.
void Shell_mult_all_imagVW_B0_FOUR
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
);


//*********************************************************************************************


/** @brief Compute \f$ (\vec{B_0} \cdot \vec{K}) \Im [\vec{A}(K) \cdot \conj{\vec{B}(K) ] \f$
 *			over a 3D ring.
 *
 * @param  N[]  The size of the array A.
 * @param  B0	The mean magnetic field.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bx  x-component of the vector field \f$ \vec{B} \f$.
 * @param  By  y-component of the vector field \f$ \vec{B} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{B} \f$.
 * @param  ring_shell_radius_array The shell radii of the rings are stored here.
 * @param  sector_angle_array The sector angles are stored here:
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  result(n,m) = \f$ (\vec{B_0} \cdot \vec{K}) 
 *							\Im [\vec{A}(K) \cdot \conj{\vec{B}(K) ] \f$  for ring(n,m).
 *
 */ 
void Local_ring_mult_all_imagVW_B0_FOUR
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> ring_shell_radius_array,  Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
);


/// Sum over local rings and send to the master node.
void Ring_mult_all_imagVW_B0_FOUR
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> ring_shell_radius_array, Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
);

//*********************************************************************************************


/** @brief Compute \f$ (\vec{B_0} \cdot \vec{K}) \Im [\vec{A}(K) \cdot \conj{\vec{B}(K) ] \f$
 *			over a 3D cylindrical ring.
 *
 * @param  N[]  The size of the array A.
 * @param  B0	The mean magnetic field.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bx  x-component of the vector field \f$ \vec{B} \f$.
 * @param  By  y-component of the vector field \f$ \vec{B} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{B} \f$.
 * @param  cylinder_shell_radius_array The shell radii of the cylindrical rings are stored here.
 * @param  cylinder_kpll_array The slab coordinates are stored here:
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  result(n) = \f$ (\vec{B_0} \cdot \vec{K}) 
 *				\Im [\vec{A}(K) \cdot \conj{\vec{B}(K) ] \f$  for cylindrical ring(n,m).
 */
void Local_cyl_ring_mult_all_imagVW_B0_FOUR
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_result, 
	DP kfactor[]
);

/// Sum over local cylindrical rings and send to the master node.
void Cyl_ring_mult_all_imagVW_B0_FOUR
(
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
);

//*********************************************************************************************


DP Local_shell_mult_vorticity_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 );

DP Shell_mult_vorticity_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 );

DP Local_shell_mult_vector_potential_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 );

DP Shell_mult_vector_potential_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 );

void Local_shell_mult_vorticity_all_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 );

void Shell_mult_vorticity_all_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 );


void Local_shell_mult_vector_potential_all_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 );

void Shell_mult_vector_potential_all_FOUR
(
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 );




#endif

//****************************  End of four_ET.h    *******************************************






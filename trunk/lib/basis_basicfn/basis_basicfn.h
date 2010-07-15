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

/*! \file	basis_basicfn.h
 * 
 * @brief basic functions declarations that are common to all basis functions.
 *
 *	shell(n) = \f$  K' \in (R^{sh}(n-1), R^{sh}(n) ] \f$.  
 *	with typical \f$ R^{sh} = 0,2,4,8,..., R_{max},\infty \f$.
 *
 *	ring(n,m) = \f$ K' = (R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ] \f$.
 *	with typical \f$ R^{sh} = 0,4,8,..., R_{max},\infty \f$ and
 *		\f$ \Theta = 0:\Theta_{max}\f$ divided in m intervals.  \f$\Theta_{max} \f$ could be
 *		\f$ \pi \f$ or \f$ \pi/2 \f$ (see function Get_max_anis_theta).
 *
 * @version 4.0 Parallel version
 * @author  M. K. Verma
 * @date	Sept 2008
 * @bug		No known bug
 */ 



#ifndef _BASIC_BASISFN_H
#define _BASIC_BASISFN_H
  
#include <blitz/array.h>
#include <complex>
#include <cmath>

#include <mpi.h>
#include <fftw3-mpi.h>

#include "../basis_basicfn/basis_basicfn_inline.h"


using namespace blitz ;


//*********************************************************************************************

/** @brief Transpose(A) -> Atr.  Only x-y components.
 *
 * Steps in each processor
 * (1) Take a xy plane Axy for a given z (in each processor).
 * (2) Transpose Axy in each processor. Axy -> Axy_tr.
 * (3) MPI_Alltoall(Axy_tr -> temp2)
 * (4) Transpose each transmitted piece in each processor.
 *
 * @param  A() complex array
 * @param  N[] size of A: (N1, N2, N3/2+1).
 * 
 * @return   Atr():  Transpose of A along xy directions.
 *
 * @sa figure ..
 *
 */
void Transpose_array(int N[], Array<complx,3> A, Array<complx,3> Atr);



//*********************************************************************************************

/** @brief Inverse_transpose(Atr) -> A.  Only x-y components.
 *
 * Steps in each processor
 * (1) Take a xy plane Atr_xy for a given z (in each processor).
 * (2) Transpose Atr_xy in each processor. Atr_xy -> Axy.
 * (3) MPI_Alltoall(Axy -> temp1)
 * (4) Transpose each transmitted piece in each processor.
 *
 * @param  Atr() complex array
 * @param  N[] size of A: (N1, N2, N3/2+1).
 * 
 * @return   A():  Transpose of Atr along xy directions.
 *
 * @sa figure ..
 *
 */
void Inverse_transpose_array(int N[], Array<complx,3> Atr, Array<complx,3> A);



//*********************************************************************************************

/** @brief Returns the index of the shell that contains \f$ \vec{K} \f$ with 
 *			\f$ |\vec{K}| \f$ = kkmag.  
 *			For flux and shell-to-shell calculations.
 *
 * @param  kkmag	magnitude of \f$ \vec{K} \f$.
 * @param  shell_radius_array	Array containing the shell radii; 
 *			for example [0,2,4,...,Rmax, infty].
 * 
 * @return   shell_index	The index of the shell that contains the vector \f$ \vec{K} \f$
 * 
 * @note	shell(n) = \f$  K' \in (R^{sh}(n-1), R^{sh}(n) ] \f$.  
 *	with typical \f$ R^{sh} = 0,2,4,8,..., R_{max},\infty \f$ <BR>
 *	shell(0) contains only K=0 mode. shell(1) contains (0:1], and so on..
 *
 */
int Get_shell_index(DP kkmag, Array<DP, 1> shell_radius_array);



//*********************************************************************************************


/** @brief Returns the index of the sector that contains \f$ \vec{K} \f$ that makes an angle
 * \f$ \theta \f$ with the anisotropy axis.
 *	For anisotropic ring spectrum.
 *
 * @param  theta	\f$ atan(K_{\perp}/K_{||}) \f$.
 * @param	sector_angle_array  Array containing the sector angles 
 *			\f$ \Theta = 0:\Theta_{max}\f$ divided in m intervals.
 *
 * @return	sector_index The index of the sector that contains the vector \f$ \vec{K} \f$.
 *
 * @note  sector(1) = [0:theta(1)]; sector(2)=(theta(1):theta(2)]; ..
 *			sector(nsector) = (theta(last-1):theta(last)).
 *
 */					
int Get_sector_index(DP theta, Array<DP, 1> sector_angle_array);


//*********************************************************************************************

/** @brief Returns the index of the ring that contains \f$ \vec{K} \f$ with 
 *				\f$ |\vec{K}| \f$ = kkmag and that makes an angle
 *				\f$ \theta \f$ with the anisotropy axis.
 *
 *	For anisotropic ring-to-ring calculation.
 *
 * @param  kkmag	magnitude of \f$ \vec{K} \f$.
 * @param  theta	\f$ atan(K_{\perp}/K_{||}) \f$.
 * @param  shell_radius_array	Array containing the shell radii; 
 *			for example [0,2,4,...,Rmax, infty].
 * @param	sector_angle_array  Array containing the sector angles = 
 *			\f$ \Theta = 0:\Theta_{max}\f$ divided in m intervals.
 * 
 * @return  shell_index The index of the shell that contains the vector \f$ \vec{K} \f$.
 * @return	sector_index The index of the sector that contains the vector \f$ \vec{K} \f$.
 * 
 * @note	shell(n) = \f$  K' \in (R^{sh}(n-1), R^{sh}(n) ] \f$.  
 *			with typical \f$ R^{sh} = 0,2,4,8,..., R_{max},\infty \f$.
 *
 * @note 	ring(n,m) = \f$ K' = (R^{ring}(n-1), \Theta(m-1); R^{ring}(n), \Theta(m) ] \f$.
 *			with typical \f$ R^{sh} = 0,4,8,..., R_{max},\infty \f$ and
 *			\f$ \Theta = 0:\Theta_{max}\f$ divided in m intervals. <BR>
 *			ring(0,1) contains only K=0  mode; other ring(0,i) contains no modes;
 *			ring(1,1) contains \f$ K \in (0,4] \f$ and \f$ \theta \in [0,\theta(1)] \f$; 
 *			ring(1,2) contains \f$ K \in (0,4] \f$ and \f$ \theta \in (\theta(1),\theta(2)] \f$
 *			and so on..
 *
 */
void Compute_ring_index
(
	DP kkmag, 
	DP theta, 
	Array<DP, 1> shell_radius_array, 
	Array<DP, 1> sector_angle_array, 
	int& shell_index, 
	int& sector_index
);
	
//*********************************************************************************************

/*! @brief Returns slab index that contains \f$ \vec{K}= (K_{||}, K_{\perp}) \f$	
 *
 * @param	kkpll	\f$ K_{||} \f$
 * @param	kkperp	\f$ K_{\perp} \f$
 * @param	cylinder_kpll_arrray Array containing kpll
 *
 * @return slab index 
 *
 * @note	slab(n) = \f$ K_{\pll} = (H(n-1), H(n)] \f$.  <BR>
 *			slab(1) contains \f$ k_{||} \in [0,H(1)] \f$ (note including  \f$  k_{||}=0 \f$.
 *			slab(2) contains \f$ k_{||} \in (H(1),H(2)] \f$
 */
int Get_slab_index(DP kkpll, DP kkperp, Array<DP,1> cylinder_kpll_array);		



//*********************************************************************************************

/*! @brief Returns slab index that contains \f$ \vec{K}= (K_{||}, K_{\perp}) \f$	
 *
 * @param	kkpll		\f$ K_{||} \f$
 * @param	kkperp		\f$ K_{\perp} \f$
 * @param	cylinder_shell_radius_array Array containing radii of the cylinders.
 * @param	cylinder_kpll_arrray Array containing kpll.
 *
 * @return	shell_index
 * @return	slab_index
 *
 * @note 	ring(n,m) = \f$ K' = (R^{ring}(n-1), H(m-1); R^{ring}(n), H(m) ] \f$. <BR>
 *			\f$ K_{||}=0 \f$ belong to slab(1).
 *			ring(0,1) contains only K=0 mode;  other ring(0,i) contains no modes;
 *			ring(1,1) contains \f$ K \in (0,4] \f$ and \f$ k_{||} \in [0,H(1)] \f$;
 *			ring(1,2) contains \f$ K \in (0,4] \f$ and \f$ k_{||} \in (H(1),H(2)] \f$
 */
void Compute_cylinder_ring_index
(
	DP kkpll, 
	DP kkperp, 
	Array<DP,1> cylinder_shell_radius_array,
	Array<DP,1> cylinder_kpll_array, 
	int& shell_index, 
	int& slab_index
);					


//*********************************************************************************************

/*! @brief computer Fourier index for the last wavenumber given real-space index.
 *
 * @param	i_last  Real-space index for last wavenumber (e.g., kz for 3D).
 *
 * @return  k_last  Wavenumber index for the last wavenumber
 * @return  is_real return 1 if real-part else 0.
 */
void Real_to_fourier_index(int i_last, int& k_last, int& is_real);

#endif

//***********************************  End of basis_basicfn.h  ********************************






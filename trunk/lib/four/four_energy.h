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




/*! \file four_energy.h 
 * 
 * @brief Contains useful functions for computing total energy, energy spectrum etc. of a 
 *			complex array in 1D, 2D, 3D.
 * 
 * Data array is complex A:  \f$ A[0:N_1/2] \f$ (1D); \f$ A[0:N_1-1, 0:N_2/2] \f$ (2D);  
 *		\f$ A[0:N_1-1, 0:N_2-1, 0:N_3/2] \f$ (3D). <BR>
 * 
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the grid index 
 *		\f$ \vec{i} \f$.
 * Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber \f$ \vec{k} \f$ 
 *		using \f$ K_i = k_i * f_i \f$
 * where \f$ f_i \f$ is the kfactor[i].  In the program \f$ K_i \f$ is written as kki.
 *
 * Modal energy  = \f$ E(k) =  |A(\vec{k})|^{2}/2 \f$   <BR>
 * Dissipation rate = \f$ D(k) = 2 K^2 E(k) \f$. <BR>
 *
 * Total energy = \f$ E = \sum_k |A(\vec{k})|^2 /2 \f$.  Since only half the array elements 
 *		are stored, we multiply each term |A(\vec{k})|^2 by a factor Multiplicity_factor_FOUR
 *		(\f$ \vec{i} \f$, N).  <BR>
 *
 * We also compute \f$ E^{AB} = \sum_k \Re(A*B^{*})/2 \f$ using the above procedure.  <BR>
 * 
 * Shell energy spectrum = \f$ Sk(s,n) = \sum[ K'^n |A(\vec{K}')|^{2} /2 ] \f$ with 
 *		\f$ s-1 < K' <= s) \f$. <BR>
 *		\f$ Sk(0,0) = |A(0)|^{2} /2 \f$ and \f$ Sk(0,n>0) = 0 \f$.
 *
 * For ring spectrum  \f$ Sk(s,m) = \sum[ K'^n E(K') ] \f$  for 
 *		\f$ \vec{K'} \in (s-1, \Theta(m-1); s, \Theta(m) ] \f$ (s>0, m>0) with 
 *		\f$ \Theta = 0:\Theta_{max} \f$ divided in m sections.  <BR>
 *		\f$ Sk(0,1) = |A(0)|^{2} /2 \f$, and \f$ Sk(0,i) = 0 \f$ for other i's. <BR>
 *		\f$ Sk(1,1) \f$ contains the contributions from \f$ K' \in (0,1] \f$ and
 *		\f$ \theta \in [0:\theta(1)] \f$; 
 *		\f$ Sk(1,2) \f$ contains the contributions from \f$ K' \in (1,2] \f$ and
 *		\f$ \theta \in (\theta(1):\theta(2)] \f$, etc. 
 *
 * Ring energy spectrum has two components: Azimuthal (toriodal) and Polar (poloidal). <BR>
 *
 * Entropy = \f$ \sum p_i log(1/p_i) \f$ where \f$ p_i = E_i/E_{total} \f$ with \f$ E_i \f$ as 
 *		the energy of the ith mode. <BR>
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 30 August 2008
 * 
 * @bug  No known bug
 */
 
 

#ifndef _FOUR_ENERGY_H
#define _FOUR_ENERGY_H

#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>

#include "four_basic.h"

using namespace blitz;

//*********************************************************************************************

/** @brief Fluctuating energy in 3D complex array A in local processor.
 *  
 * Remove the mean energy if master processor.
 * Multiply by 2 to each mode since complex conj are not stored.  
 * The modes on the xy-plane are not doubled since their complex conj are stored. <BR>
 * we count ky=-N[2]/2 modes twice since ky=-N[2]/2  modes are not stored. <BR>
 * kx=-N[1]/2 & ky=-N[2]/2  modes are not stored. So 
 * If numproc == 1 (sequential), we count modes with l1=N[1]/2 (kx=N[1]/2) twice.
 * Else, for my_id=numprocs/2, double modes with l1=0 (kx=N[1]/2).
 *
 * @param  N[]  The size of the array A.
 * @param  A	Local complex array.
 * 
 * @return local_E = \f$ \sum_k |A(\vec{k})|^2 /2 \f$.
 */
DP Get_local_energy_FOUR(string alias_switch, int N[], Array<complx,3> A);

/// Returns sum all local_E's.
DP Get_total_energy_FOUR(string alias_switch, int N[], Array<complx,3> A);

//*********************************************************************************************


/** @brief Fluctuating energy in 3D complex array A in local processor.
 *  
 * Remove the mean energy if master processor.
 * Multiply by 2 to each mode since complex conj are not stored.  
 * The modes on the xy-plane are not doubled since their complex conj are stored. <BR>
 * we count ky=-N[2]/2 modes twice since ky=-N[2]/2  modes are not stored. <BR>
 * kx=-N[1]/2 & ky=-N[2]/2  modes are not stored. So 
 * If numproc == 1 (sequential), we count modes with l1=N[1]/2 (kx=N[1]/2) twice.
 * Else, for my_id=numprocs/2, double modes with l1=0 (kx=N[1]/2).
 *
 * @param  N[]  The size of the arrays A & B.
 * @param  A	Local complex array.
 * @param  B	Local complex array.
 * 
 * @return local_E = \f$ \Real(A*conj(B))/2 \f$.
 */
DP Get_local_energy_FOUR(string alias_switch, int N[],  Array<complx,3> A, Array<complx,3> B);

/// Returns sum all local_E's.
DP Get_total_energy_FOUR(string alias_switch, int N[],  Array<complx,3> A, Array<complx,3> B);

//*********************************************************************************************

/** @brief Computes \f$ \sum  K^n |A(k)|^2/2 \f$ in local processor.
 *  
 * @param  N[]  The size of the array A.
 * @param  A  Local complex array.
 * @param  n  exponent 
 * 
 * @return local_Sn = \f$ \sum  K^n |A(k)|^2/2 \f$.
 *
 * @sa Get_local_energy_FOUR(int N[], Array<complx,3> A)
 */
DP Get_local_Sn_FOUR(string alias_switch, int N[], Array<complx,3> A, DP n, DP kfactor[]);

/// Returns sum all local_Sn's.
DP Get_total_Sn_FOUR(string alias_switch, int N[], Array<complx,3> A, DP n, DP kfactor[]);


//*********************************************************************************************

/** @brief local shell spectrum \f$ Sk(s,n) = \sum[ K'^n |A(\vec{K}')|^{2} /2 ]  \f$  in the local processor.
 *  
 * Multiply by 2 to each mode since complex conj are not stored.  
 * The modes on the xy-plane are not doubled since their complex conj are stored.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.   
 * @param  n  exponent 
 * 
 * @return  local_Sk(s,n) = \f$  \sum[ K'^n |A(\vec{K}')|^{2} /2 ] \f$ with \f$ s <= K' <(s+1) \f$.
 * @return  local_Sk_count(s,n) number of wavenumber lattice points in the shell \f$ s <= K' <(s+1) \f$.
 */
void Compute_local_shell_spectrum_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, 
	DP n, 
	Array<DP,1> local_Sk, 
	Array<DP,1> local_Sk_count,  
	DP kfactor[]
);



/** @brief  Compute the total shell spectrum by summing all the local shell spectra.
 *  
 * The shells near the edge are semifilled, so the spectrum for these shells are extrapolated.
 * The total number of lattice points on these shells are sum of local_Sk_count.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.   
 * @param  n  exponent 
 * 
 * @return  Sk(s,n) = sum of local_Sk.  
 * The edge shells are appropriately extrapolated.
 * 
 * @sa Compute_local_shell_spectrum_FOUR(int N[],  Array<complx,3> A, DP n, Array<DP,1> local_Sk, Array<DP,1> local_Sk_count,  DP kfactor[])
 */
void Compute_shell_spectrum_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, 
	DP n, 
	Array<DP,1> Sk,  
	DP kfactor[]
);



/** @brief local shell spectrum real(A(k).B(k)*)/2  in the local processor.
 *  
 * Multiply by 2 to each mode since complex conj are not stored.  
 * The modes on the xy-plane are not doubled since their complex conj are stored.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array. 
 * @param  B  a complex array.  
 * @param  n  exponent 
 * 
 * @return  local_Sk(s,n) = \f$  \sum[ K'^n \Re(A(\vec{K}') * B^{*}(\vec{K}')) /2] \f$ with \f$ s <= K' <(s+1) \f$.
 * @return  local_Sk_count(s,n) number of wavenumber lattice points in the shell \f$ s <= K' <(s+1) \f$.
 */								
void Compute_local_shell_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> A,  Array<complx,3> B,
	DP n, 
	Array<DP,1> local_Sk, 
	Array<DP,1> local_Sk_count, 
	DP kfactor[]
);


/** @brief  Compute the total shell spectrum by summing all the local shell spectra.
 *  
 * @return  Sk(s,n) = sum of local_Sk = \f$  \sum[ K'^n \Re(A(\vec{K}') * B^{*}(\vec{K}')) /2] \f$ with \f$ s <= K' <(s+1) \f$.  
 * The edge shells are appropriately extrapolated.
 * 
 * @sa Compute_local_shell_spectrum_FOUR(int N[],  Array<complx,3> A, Array<complx,3> B, DP n, 
 *				Array<DP,1> local_Sk, Array<DP,1> local_Sk_count,  DP kfactor[])
 */
void Compute_shell_spectrum_FOUR
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP n, 
	Array<DP,1> Sk,  
	DP kfactor[]
);							



//*********************************************************************************************


/** @brief Computes local sum of helicity1, helicity2, dissipation_helicity1, dissipation_helicity2.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  local_helicity1 = \f$ \sum \vec{K} . (\vec{Vr} \times \vec{Vi}) \f$; 
 * @return  local_helicity2 = \f$ \sum [\vec{K} . (\vec{Vr} \times \vec{Vi})] /K^2 \f$.
 * @return  local_dissipation_H1 = \f$ 2*\sum [K^2 \vec{K} . (\vec{Vr} \times \vec{Vi})] \f$.
 * @return  local_dissipation_H2 = \f$ 2*sum[ \vec{K} . (\vec{Vr} \times \vec{Vi}) ] = 2*H2 \f$.
 */							
void Compute_local_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &local_helicity1, DP &local_helicity2, 
	DP &local_dissipation_H1, DP &local_dissipation_H2,
	DP kfactor[]
);


/** @brief Computes local sum of helicity1, helicity2, dissipation_helicity1, dissipation_helicity2.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  total_helicity1 = \f$ \sum \vec{K} . (\vec{Vr} \times \vec{Vi}) \f$; 
 * @return  total_helicity2 = \f$ \sum [\vec{K} . (\vec{Vr} \times \vec{Vi})] /K^2 \f$.
 * @return  total_dissipation_H1 = \f$ 2*\sum [K^2 \vec{K} . (\vec{Vr} \times \vec{Vi})] \f$.
 * @return  total_dissipation_H2 = \f$ 2*sum[ \vec{K} . (\vec{Vr} \times \vec{Vi}) ] = 2*H2 \f$.
 *				by summing local_helicity1, local_helicity2, 
 *
 * @sa Compute_local_helicity_FOUR(int N[], Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, DP kfactor[],
 *								DP &local_helicity1, DP &local_helicity2, total_dissipation_H1, and total_dissipation_H2.
 *								DP &local_dissipation_H1, DP &local_dissipation_H2)
 */	
void Compute_total_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &total_dissipation_H1, DP &total_dissipation_H2,
	DP kfactor[]
);


/** @brief Computes local shell spectrum \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return local_H1(K) = \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$ for  \f$ K \le K' K+1 \f$ (K>0).
 * @return local_H1k_count(s,n) number of wavenumber lattice points in the shell \f$ s <= K' <(s+1) \f$.
 */								
void Compute_local_shell_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3, 
	Array<DP,1> local_H1k_count, 
	DP kfactor[]
);											 


/** @brief	Computes shell spectrum \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$.
 *
 * @return	H1(K) = \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$  by summing the local sums.
 *
 * @sa	Compute_local_shell_spectrum_helicity_FOUR(int N[], Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
 *											 Array<DP,1> local_H1k, Array<DP,1> local_H1k_count, DP kfactor[])
 */	
void Compute_shell_spectrum_helicity_FOUR
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3, 
	DP kfactor[]
);


//*********************************************************************************************

/** @brief Computes local entropy of a 3D vector field A.
 *
 *  Probability is defined using energy of the mode \f$ p(\vec{k}) =  E(\vec{k})/E_{total) \f$ 
 *  with \f$ E(\vec{k}) = |\vec{A}(\vec{k})|^2/2 \f$.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 *
 * @return  local_entropy = \f$ \sum p(\vec{k}) log(1/p(\vec{k})) \f$.	
 */																																													
DP Get_local_entropy_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
);

/** @brief Computes total  entropy of a 3D vector field A by summing up the local_entropy.
 *
 * @return  entropy = \f$ \sum p(\vec{k}) log(1/p(\vec{k})) \f$ by summing up local_entropy.	
 *
 * @sa		Get_local_entropy_FOUR(int N[], Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az)
 */
DP Get_entropy_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
);


/** @brief Computes local entropy of a 3D ccalar field A.
 *
 *  Probability is defined using energy of the mode \f$ p(\vec{k}) =  E(\vec{k})/E_{total) \f$ with \f$ E(\vec{k}) = |A(\vec{k})|^2/2 \f$.
 *
 * @param  N[]  The size of the array A.
 * @param  A  3D scalar field.
 *
 * @return local_entropy = \f$ \sum p(\vec{k}) log(1/p(\vec{k})) \f$.	
 */
DP Get_local_entropy_scalar_FOUR(string alias_switch, int N[], Array<complx,3> A);

/** @brief Computes total entropy of a 3D ccalar field A by summing local_entropy.
 *
 * @return total_entropy = \f$ \sum p(\vec{k}) log(1/p(\vec{k})) \f$ by summing up local_entropy.
 *
 * @sa		Get_local_entropy_scalar_FOUR(int N[], Array<complx,3> A)
 */
DP Get_entropy_scalar_FOUR(string alias_switch, int N[], Array<complx,3> A);	


//*********************************************************************************************


/** @brief Computes local ring spectrum  of a 3D vector fields A.B.
 *
 * Craya decomposition is used. <BR>
 * \f$ \vec{e}_1 = \vec{K} \times \vec{n} /|\vec{K} \times \vec{n}| \f$ <BR>
 * \f$ \vec{e}_2 = \vec{K} \times \vec{e}_1 /K \f$ <BR>
 * where \f$ \vec{n} \f$ unit vector along the anisotropy direction.
 * 
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  sector_angle_array The sector angles are stored here:\f$ \Theta = 0:\Theta_{max}\f$.
 *
 * @return  \f$ local_S1k(s,m) = \sum[ K'^n |\vec{A}(\vec{K'}) \cdot \vec{e}_1(K')|^2 /2 ] \f$ 
 * @return  \f$ local_S2k(s,m) = \sum[ K'^n |\vec{A}(\vec{K'}) \cdot \vec{e}_2(K')|^2 /2 ] \f$
 *				for \f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).
 */
void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
);


/// Total ring spectrum by summing the local ring spectrum.
void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);						



/** @brief Computes local ring spectrum  of a 3D vector fields A.B.
 *
 * Craya decomposition is used. <BR>
 * \f$ \vec{e}_1 = \vec{K} \times \vec{n} /|\vec{K} \times \vec{n}| \f$ <BR>
 * \f$ \vec{e}_2 = \vec{K} \times \vec{e}_1 /K \f$ <BR>
 * where \f$ \vec{n} \f$ unit vector along the anisotropy direction.
 * 
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  sector_angle_array The sector angles are stored here:\f$ \Theta = 0:\Theta_{max}\f$.
 *
 * @return  \f$ local_S1k(s,m) = \sum[ K'^n |\vec{A}(\vec{K'}) \cdot \vec{e}_1(K')|^2 /2 ] \f$ 
 * @return  \f$ local_S2k(s,m) = \sum[ K'^n |\vec{A}(\vec{K'}) \cdot \vec{e}_2(K')|^2 /2 ] \f$
 *				for \f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).
 */
void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
);


/// Total ring spectrum by summing the local ring spectrum.
void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);


/** @brief Compute local ring spectrum of 3D scalar field A.
 *
 * @param  N[]  The size of the array A.
 * @param  F   complex scalar field.  
 * @param  sector_angle_array The sector angles are stored here: \f$\Theta = 0:\Theta_{max}\f$.
 *
 * @return  \f$ local_Sk(s,m) = \sum[ K'^n |F(\vec{K'})|^2   /2 ] \f$  
 *				for \f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).		
 */
void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
);


/** @brief Compute local ring spectrum of 3D scalar field A.
 *
 * @return  \f$ Sk(s,m) = \sum[ K'^n |F(\vec{K'})|^2   /2 ] \f$  
 *			by summing up local ring spectrum 	local_Sk(s,m)
 */
void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
);


/** @brief Compute local ring spectrum F*G of 3D scalar field F & G.
 *
 * @param  N[]  The size of the array A.
 * @param  F,G   complex scalar fields.  
 * @param  sector_angle_array The sector angles are stored here: \f$\Theta = 0:\Theta_{max}\f$.
 *
 * @return  \f$ local_Sk(s,m) = \sum[ K'^n \Re(F(\vec{K'}) * G^*(\vec{K'}))   /2 ] \f$  
 *				for \f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).		
 */
void Compute_local_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
);


/** @brief Compute local ring spectrum F*G of 3D scalar field F & G.
 *
 * @return  \f$ Sk(s,m) = \sum[ K'^n \Re(F(\vec{K'}) * G^*(\vec{K'}))   /2 ] \f$    
 *			by summing up local ring spectrum 	local_Sk(s,m)
 */
void Compute_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> Sk, 
	DP kfactor[]
);



/** @brief Computes local ring helicity spectrum  of a 3D vector field A.
 * 
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  sector_angle_array The sector angles are stored here:\f$ \Theta = 0:\Theta_{max}\f$.
 *
 * @return  \f$ local_H1k(s,m) = \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$ 
 *				for \f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).
 */	 									
void Compute_local_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array,
	Array<DP,2> local_H1k, 
	DP kfactor[]
);	


/** @brief Computes local ring helicity spectrum  of a 3D vector field A.
 *
 * @return  \f$ H1k(s,m) = \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$ 
 *				by summing local_H1k.
 */	
void Compute_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array,
	Array<DP,2> H1k, 
	DP kfactor[]
);



//*********************************************************************************************

// Cylinders
																	

/** @brief Computes local cylinder_ring spectrum  of a 3D vector field A.
 *
 * \f$ \vec{e}_1 : \f$ Along the anisotropy direction <BR>
 * \f$ \vec{e}_2 : \f$ Perpendicular to the anisotropy direction.
 * 
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$. 
 * @param  cylinder_kpll_array: [-min_kpll, ..., max_kpll]. 
 *
 * @return  \f$ S1k(s,m) = \sum |\vec{A}(\vec{K'}) \cdot \vec{e}_1(K')|^2 /2 \f$.
 * @return  \f$ S2k(s,m) = \sum |\vec{A}(\vec{K'}) \cdot \vec{e}_2(K')|^2 /2 \f$
 *				for \f$ K' = (s-1, H(m-1); s, H(m) ] \f$ (s>0, m>0).
 */
void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
);


/// Compute total cylindrical ring spectrum by summing up local cylindrical ring spectrum
/// for vector field A.
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);




/** @brief Computes local cylinder_ring spectrum  of 3D vector fields A.B.
 *
 * \f$ \vec{e}_1 : \f$ Along the anisotropy direction <BR>
 * \f$ \vec{e}_2 : \f$ Perpendicular to the anisotropy direction.
 * 
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 *
 * @param  Bx  x-component of vector field \f$ \vec{A} \f$.
 * @param  By  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{A} \f$.
 * @param  cylinder_kpll_array: [-min_kpll, ..., max_kpll].
 *
 * @return  \f$ S1k(s,m) = \sum (\vec{A}(\vec{K'}) \cdot \vec{e}_1(K')) 
						(\vec{B}^{*}(\vec{K'}) \cdot \vec{e}_1(K')) /2 \f$ 
 *				for \f$ K' = (s-1, H(m-1); s, H(m) ] \f$  (s>0, m>0).
 * @return  \f$ S2k(s,m) = \sum (\vec{A}(\vec{K'}) \cdot \vec{e}_2(K')) 
						(\vec{B}^{*}(\vec{K'}) \cdot \vec{e}_2(K')) /2 \f$.		
 */		
void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
);


/// Compute total cylindrical ring spectrum by summing up local cylindrical ring spectrum
/// for vector fields A.B
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);


/** @brief Compute  cylinder_ring spectrum of 3D scalar field A.
 *
 * @param  N[]  The size of the array A.
 * @param  F   complex scalar field.
 * @param cylinder_kpll_array: [-min_kpll, ..., max_kpll].
 *
 * @return  \f$ Sk(s,m) = \sum |F(\vec{K'})|^2 /2 \f$  
 *			for \f$ K' = (s-1, H(m-1); s, H(m) ] \f$  (s>0, m>0).		
 */	
void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> local_Sk, 
	DP kfactor[]
);


/// Compute total cylindrical ring spectrum by summing up local cylindrical ring spectrum
/// for scalar field F
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> Sk, 
	DP kfactor[]
);


/** @brief Compute local cylinder_ring spectrum F.G for 3D scalar fields F & G.
 *
 * @param  N[]  The size of the array F.
 * @param  F,G   complex scalar field.
 * @param cylinder_kpll_array: [-min_kpll, ..., max_kpll].
 *
 * @return  \f$ Sk(s,m) = \sum |F(\vec{K'})|^2 /2 \f$  
 *			for \f$ K' = (s-1, H(m-1); s, H(m) ] \f$  (s>0, m>0).		
 */	
void Compute_local_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G,
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> local_Sk, 
	DP kfactor[]
);


/// Compute total cylindrical ring spectrum by summing up local cylindrical ring spectrum
/// for scalar fields F & G.
void Compute_cylinder_ring_spectrum_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<complx,3> G,
	Array<DP, 1> cylinder_kpll_array,
	DP n,   
	Array<DP,2> Sk, 
	DP kfactor[]
);



/** @brief Computes local cylinder_ring helicity spectrum  for vector field A.
 * 
 * @param  N[]  The size of the array A.
 * @param  A  complex vector field.
 * @param  cylinder_kpll_array: [-min_kpll, ..., max_kpll].
 *
 * @return  \f$ H1k(s,m) = \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$ 
 *				for \f$ K' = (s-1, H(m-1); s, H(m) ] \f$ (s>0, m>0).	
 */
void Compute_local_cylinder_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_H1k, 
	DP kfactor[]
);


/// Compute total helicity ring spectrum for cylindrical geometry by summing up 
/// local cylindrical ring spectrum for vector field A.
void Compute_cylinder_ring_spectrum_helicity_FOUR
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> H1k, 
	DP kfactor[]
);



//*********************************************************************************************

/** @brief Compute local spectrum imag(\f$ (\vec{B0} \cdot \vec{K'}) \vec{A}.\vec{B}^{*} \f$) 
 *					over all 2D shells.
 *
 * @param  N[]  The size of the array A.
 * @param  B0  Mean magnetic field.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 *
 * @param  Bx  x-component of vector field \f$ \vec{A} \f$.
 * @param  By  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Sk  shell-spectrum for B0 energy transfer
 *
 * @return  Sk(n) =\f$  \sum[ (\vec{B0} \cdot \vec{K'}) 
 *						\Im(\vec{A}(\vec{K}') * \vec{B}^{*}(\vec{K}')) ] \f$ for 
 *		  \f$ s-1 < K' <= s \f$.
 */ 
void Compute_local_imag_shell_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> local_Sk, 
	DP kfactor[]
);



/// Compute total spectrum by summing up local spectrum.
void Compute_imag_shell_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP,1> Sk, 
	DP kfactor[]
);



/** @brief Compute spectrum real(\f$ (\vec{B0} \cdot \vec{K'}) \vec{A}.\vec{B}^{*} \f$) 
 *					over all 3D ring.
 *
 * @param  N[]  The size of the array A.
 * @param  B0  Mean magnetic field.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 *
 * @param  Bx  x-component of vector field \f$ \vec{A} \f$.
 * @param  By  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Sk  ring-spectrum for B0 energy transfer
 *
 * @return  Sk(n) =\f$  \sum[ (\vec{B0} \cdot \vec{K'}) 
 *						\Im(\vec{A}(\vec{K}') * \vec{B}^{*}(\vec{K}')) ] \f$ for 
 *				for \f$ K' = (s-1, \Theta(m-1); s, \Theta(m) ] \f$  (s>0, m>0).
 */
void Compute_local_imag_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
);


/// Compute total spectrum by summing up local spectrum.
void Compute_imag_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> sector_angle_array, 
	Array<DP,2> Sk, 
	DP kfactor[]
);

/** @brief Compute local spectrum real(\f$ (\vec{B0} \cdot \vec{K'}) \vec{A}.\vec{B}^{*} \f$) 
 *					over all 3D cylindrical rings.
 *
 * @param  N[]  The size of the array A.
 * @param  B0  Mean magnetic field.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 *
 * @param  Bx  x-component of vector field \f$ \vec{A} \f$.
 * @param  By  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Bz  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Sk cylindrical-ring-spectrum for B0 energy transfer
 *
 * @return  Sk(n) =\f$  \sum[ (\vec{B0} \cdot \vec{K'}) 
 *						\Im(\vec{A}(\vec{K}') * \vec{B}^{*}(\vec{K}')) ] \f$ for 
 *				\f$ K' = (R^{ring}(n-1), H(m-1); R^{ring}(n), H(m) ] \f$.
 */
void Compute_local_imag_cylinder_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
);


/// Compute total spectrum by summing up local spectrum.
void Compute_imag_cylinder_ring_spectrum_B0_FOUR
(
	string alias_switch,
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> Sk, 
	DP kfactor[]
);

							
																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																											
#endif

//*****************************  End of four_energy.h *****************************************	








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

/*! \file scft_energy.h 
 * 
 * @brief Contains useful functions for computing total energy, energy spectrum etc. 
 *			of a Sin/CosFOUR transformed array.
 * 
 * Data array is real \f$ A[0:N_1/2] \f$ (1D). <BR>
 * Data array is complex A:  \f$ A[0:N_1-1, 0:N_2/2] \f$ (2D);  
 *			\f$ A[0:N_1-1, 0:N_2-1, 0:N_3/2] \f$ (3D). <BR>
 * 
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the grid index 
 *		\f$ \vec{i} \f$. Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber 
 *		\f$ \vec{k} \f$ using \f$ K_i = k_i * f_i \f$ where \f$ f_i \f$ is the kfactor[i].  
 *		In the program \f$ K_i \f$ is written as kki.
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber or 
 *		gridwavenumber depending on the switch WAVENOACTUAL or WAVENOGRID.   
 *		Actual wavenumber is preferred.
 *
 * Total energy  = \f$ E =  |A(0, \vec{k}_{\perp})|^{2}/2  
 *							+ \sum_{k_x>0}  |A(\vec{k})|^2  \f$.  <BR>
 *		Note that kx=0 Sin modes are zero. <BR>
 *		Since only half the array elements are stored, we multiply each term 
 *		|A(\vec{k})|^2 by a factor Multiplicity_factor_SCFT(\f$ \vec{i} \f$, N).  <BR>
 *
 * We also compute \f$ E^{AB} \propto \Re(A*B^{*}) \f$ using the above procedure.  <BR>
 *
 * Modal energy  = \f$ E(k) =  |A(k_x, \vec{k}_{\perp})|^{2} \f$  if  \f$ (k_x > 0) \f$.  <BR>
 * Modal energy  = \f$ E(k) =  |A(0, \vec{k}_{\perp})|^{2}/2 \f$  if  \f$ (k_x = 0) \f$.  <BR>
 * Dissipation rate = \f$ D(k) = 2 K^2 E(k) \f$. <BR>
 * Note that kx=0 Sin modes are zero. <BR>
 * 
 * Shell energy spectrum = \f$ Sk(s,n) = \sum[ K'^n E(K') ] \f$ with 
 *							\f$ s-1 < K' <= s \f$. <BR>
 *
 * For ring spectrum  \f$ Sk(s,m) = \sum[ K'^n E(K') ] \f$  for \f$ K' 
 *		= [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0) with \f$ s-1 < K' <= s \f$ and 
 *		\f$ \Theta = [0:\Theta_{max}] \f$ divided into m sectors.
 *
 * Ring energy spectrum has two components: Azimuthal (toriodal) and Polar (poloidal).  <BR>
 *
 * Entropy = = \f$ \sum p(k) log(1/p(k)) \f$ where probability  \f$ p(k) =  E(k)/E_{total) \f$ 
 *		with \f$ E(k) = |A(k)|^2/2 \f$.
 *
 * @sa scft_energy.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 30 August 2008
 * 
 * @bugs  No known bug
 */

#ifndef _SCFT_ENERGY_H
#define _SCFT_ENERGY_H

#include <blitz/array.h>
#include <complex>
#include <cmath>

#include <mpi.h>
#include <fftw3-mpi.h>

#include "scft_basic.h"


using namespace blitz;

//*********************************************************************************************

/** @brief Fluctuating energy of a complex array A for ky = ny plane.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * 
 * @return  \f$ E =  |A(0, \vec{k}_{\perp})|^{2}/2  + \sum_{k_x>0}  |A(\vec{k})|^2  \f$.
 */
DP Get_local_energy_XZ_plane_SCFT(string alias_switch, int N[], Array<complx,3> A, int ny);

/** @brief Fluctuating energy in 3D complex array A in local processor.
 *  
 * Computed using Get_local_energy_XZ_plane_SCFT.
 *
 * @param  N[]  The size of the array A.
 * @param  A	Local complex array.
 * 
 * @return local_E = \f$ \sum_k |A(\vec{k})|^2 /2 \f$.
 *
 * @sa Get_local_energy_XZ_plane_SCFT(int N[], Array<complx,3> A, int ny)
 */
DP Get_local_energy_SCFT(string alias_switch, int N[], Array<complx,3> A);

/// Returns sum all local_E's.
/// Remove the mean energy if master processor.
DP Get_total_energy_SCFT(string alias_switch, int N[], Array<complx,3> A);


//*********************************************************************************************

/** @brief Fluctuating energy of a complex arrays A & B for ky = ny plane.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.
 * 
 * @return  \f$ E =  \sum \Re(A(0, \vec{k}_{\perp})*conj(B(0, \vec{k}_{\perp})))/2  
		+ \sum_{k_x>0} \Re(A(\vec{k})*conj( B(\vec{k})))  \f$.
 */
DP Get_local_energy_XZ_plane_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	int ny
);

/** @brief Fluctuating energy in 3D complex array A in local processor.
 *  
 * Computed using Get_local_energy_XZ_plane_SCFT.
 *
 * @param  N[]  The size of the arrays A & B.
 * @param  A	Local complex array.
 * @param  B	Local complex array.
 * 
 * @return local_E = \f$ \Real(A*conj(B))/2 \f$.
 *
 * @sa Get_local_energy_XZ_plane_SCFT(int N[], Array<complx,3> A, Array<complx,3> B, int ny)
 */
DP Get_local_energy_SCFT(int N[],  Array<complx,3> A, Array<complx,3> B);

/// Returns sum all local_E's.
/// Remove the mean energy if master processor.
DP Get_total_energy_SCFT(string alias_switch, int N[],  Array<complx,3> A, Array<complx,3> B);


//*********************************************************************************************


/** @brief Computes \f$ \sum  K^n E(k) \f$ in local processor.
 *  
 * @param  N[]  The size of the array A.
 * @param  A  Local complex array.
 * @param  n  exponent 
 * 
 * @return local_Sn = \f$ \sum  K^n E(k) \f$.
 *
 * @sa Get_local_energy_SCFT(int N[], Array<complx,3> A)
 */
DP Get_local_Sn_SCFT(string alias_switch, int N[], Array<complx,3> A, DP n, DP kfactor[]);

/// Returns sum all local_Sn's.
DP Get_total_Sn_SCFT(string alias_switch, int N[], Array<complx,3> A, DP n, DP kfactor[]);
				

//*********************************************************************************************

/** @brief local shell spectrum \f$ Sk(s,n) = \sum[ K'^n E(K') ]  \f$  in the local processor.
 *  
 * Multiply by 2 to each mode since complex conj are not stored.  
 * The modes on the xy-plane are not doubled since their complex conj are stored.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array.   
 * @param  n  exponent 
 * 
 * @return  local_Sk(s,n) = \f$  \sum[ K'^n E(K') ] \f$ with \f$ s <= K' <(s+1) \f$.
 * @return  local_Sk_count(s,n) number of wavenumber lattice points in the shell 
 *			\f$ s <= K' <(s+1) \f$.
 */
void Compute_local_shell_spectrum_SCFT
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
 * @sa Compute_local_shell_spectrum_SCFT(int N[],  Array<complx,3> A, DP n, 
 *				Array<DP,1> local_Sk, Array<DP,1> local_Sk_count,  DP kfactor[])
 */								
void Compute_shell_spectrum_SCFT
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, 
	DP n, 
	Array<DP,1> Sk,  
	DP kfactor[]
);


/** @brief local shell spectrum \f$ K'^n E^{AB}(K') \f$  in the local processor.
 *  
 * Multiply by 2 to each mode since complex conj are not stored.  
 * The modes on the xy-plane are not doubled since their complex conj are stored.
 *
 * @param  N[]  The size of the array A.
 * @param  A  a complex array. 
 * @param  B  a complex array.  
 * @param  n  exponent 
 * 
 * @return  local_Sk(s,n) = \f$  \sum[ K'^n E^{AB}(K') ] \f$ with \f$ s <= K' <(s+1) \f$.
 * @return  local_Sk_count(s,n) number of wavenumber lattice points in the shell 
 *						\f$ s <= K' <(s+1) \f$.
 */										
void Compute_local_shell_spectrum_SCFT
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
 * @return  Sk(s,n) = sum of local_Sk = \f$  \sum[ K'^n E^{AB}(K') ] \f$ 
 *						with \f$ s <= K' <(s+1) \f$.  
 * The edge shells are appropriately extrapolated.
 * 
 * @sa Compute_local_shell_spectrum_SCFT(int N[],  Array<complx,3> A, Array<complx,3> B, DP n, 
 *				Array<DP,1> local_Sk, Array<DP,1> local_Sk_count,  DP kfactor[])
 */									
void Compute_shell_spectrum_SCFT
(
	string alias_switch, 
	int N[],  
	Array<complx,3> A, Array<complx,3> B, 
	DP n, 
	Array<DP,1> Sk,  
	DP kfactor[]
);								

//*********************************************************************************************

/** @brief local  helicity1, helicity2, dissipation_helicity1, dissipation_helicity2  
 *			in the local processor.
 *  
 * Given SCFT vector \f$ (V_1 , V_2, V_3) \f$, the Fourier space vector is 
 *			\f$ (V_1/i , V_2, V_3) \f$.
 * Note that the first component of the velocity vector SinFOUR transformed.  
 * We are using RB free-slip basis functions as standard for this calculation. 
 *
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return  total_helicity1 = \f$ \sum \vec{K} . (\vec{Vr} \times \vec{Vi}) \f$; 
 * @return  total_helicity2 = \f$ \sum [\vec{K} . (\vec{Vr} \times \vec{Vi})] /K^2 \f$.
 * @return  total_dissipation_H1 = \f$ 2*\sum [K^2 \vec{K} . (\vec{Vr} \times \vec{Vi})] \f$.
 * @return  total_dissipation_H2 =\f$ 2*sum[ \vec{K} . (\vec{Vr} \times \vec{Vi}) ] = 2*H2 \f$.
 */
void Compute_local_helicity_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &local_helicity1, DP &local_helicity2, 
	DP &local_dissipation_H1, DP &local_dissipation_H2,
	DP kfactor[]
);	
								
/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.																									
void Compute_total_helicity_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP &total_helicity1, DP &total_helicity2, 
	DP &total_dissipation_H1, DP &total_dissipation_H2,
	DP kfactor[]
);


/** @brief Computes local shell spectrum 
 *         \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$ in local processor.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return H1(K) = \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) \f$ 
 *					for  \f$ K \le K' K+1 \f$ (K>0).
 * 
 * @sa Get_Modal_helicity_SCFT(int i1, int i2, int i3, int N[], Array<complx,3> Ax, 
 *				Array<complx,3> Ay, Array<complx,3> Az, DP kfactor[])
 */									
void Compute_local_shell_spectrum_helicity_SCFT
(
	string alias_switch,
	int N[], Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> local_H1k1, Array<DP,1> local_H1k2, Array<DP,1> local_H1k3, 
	Array<DP,1> local_H1k_count, 
	DP kfactor[]
);	

/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.											 										 
void Compute_shell_spectrum_helicity_SCFT
(
	string alias_switch, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP,1> H1k1, Array<DP,1> H1k2, Array<DP,1> H1k3, 
	DP kfactor[]
);

//*********************************************************************************************


/** @brief Computes local entropy of a 3D vector field A in local processor.
 *
 *  Probability is defined using energy of the mode \f$ p(\vec{k}) =  E(\vec{k})/E_{total) \f$ 
 *  with \f$ E(\vec{k}) = |\vec{A}(\vec{k})|^2/2 \f$.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 *
 * @return  entropy = \f$ \sum p(\vec{k}) log(1/p(\vec{k})) \f$.	
 */												 
DP Get_local_entropy_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
);

/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.
DP Get_entropy_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az
);

/** @brief Computes entropy of a 3D scalar field A in local processor.
 *
 *  Probability is defined using energy of the mode \f$ p(\vec{k}) 
 *			=  E(\vec{k})/E_{total) \f$ with \f$ E(\vec{k}) = |A(\vec{k})|^2/2 \f$.
 *
 * @param  N[]  The size of the array A.
 * @param  A  3D scalar field.
 *
 * @return entropy = \f$ \sum p(\vec{k}) log(1/p(\vec{k})) \f$.	
 */	
DP Get_local_entropy_scalar_SCFT(string alias_switch, int N[], Array<complx,3> A);


/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.
DP Get_entropy_scalar_SCFT(string alias_switch, int N[], Array<complx,3> A);																					 											 


/** @brief Computes local ring spectrum  of a 3D vector field A in local processor.
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
 * @return  \f$ S1k(s,m) = \sum[ K'^n |\vec{A}(\vec{K'}) \cdot \vec{e}_1(K')|^2 ] \f$.
 * @return  \f$ S2k(s,m) = \sum[ K'^n |\vec{A}(\vec{K'}) \cdot \vec{e}_2(K')|^2 ] \f$
 *				[exception for kx=0, factor =1/2] for 
 *				\f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).
 */	
void Compute_local_ring_spectrum_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_S1k, Array<DP,2> local_S2k, 
	DP kfactor[]
);
								
/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.								
void Compute_ring_spectrum_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> S1k, Array<DP,2> S2k, 
	DP kfactor[]
);	
								
/** @brief Computes local ring spectrum  of a 3D vector field A in local processor.
 *
 * Craya decomposition is used. <BR>
 * \f$ \vec{e}_1 = \vec{K} \times \vec{n} /|\vec{K} \times \vec{n}| \f$ <BR>
 * \f$ \vec{e}_2 = \vec{K} \times \vec{e}_1 /K \f$ <BR>
 * where \f$ \vec{n} \f$ unit vector along the anisotropy direction. <BR>
 *
 * Given SCFT vector \f$ (V_1 , V_2) \f$, the Fourier space vector is \f$ (V_1/i , V_2) \f$.
 * Note that the first component of the velocity vector SinFOUR transformed. 
 * We are using RB free-slip basis functions as standard for this calculation. <BR>
 * 
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  Bi  i-component of vector field \f$ \vec{B} \f$.
 * @param  sector_angle_array The sector angles are stored here:\f$ \Theta = 0:\Theta_{max}\f$. 
 *
 * @return  \f$ S1k(s,m) = \sum[ K'^n (\vec{A}(\vec{K'}) 
 *					\cdot \vec{e}_1(K')) (\vec{B}^{*}(\vec{K'}) \cdot \vec{e}_1(K')) ] \f$ 
 *					for \f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).
 *
 * @return  \f$ S2k(s,m) = \sum[ K'^n (\vec{A}(\vec{K'}) 
 *					\cdot \vec{e}_2(K')) (\vec{B}^{*}(\vec{K'}) \cdot \vec{e}_2(K'))] \f$ 
 *					[exception for kx=0, factor =1/2] .		
 */																			
void Compute_local_ring_spectrum_SCFT
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
								
/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.									
void Compute_ring_spectrum_SCFT
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
								
/** @brief Compute local ring spectrum of 3D scalar field A in local processor.
 *
 * @param  N[]  The size of the array A.
 * @param  F   complex scalar field.
 * @param  sector_angle_array The sector angles are stored here:\f$ \Theta = 0:\Theta_{max}\f$. 
 *
 * @return  \f$ Sk(s,m) = \sum[ K'^n E(K') ] \f$  for 
 *				\f$ K' = [s, \Theta(m-1); s+1, \Theta(m) ) \f$ (s>0, m>0).		
 */									
void Compute_local_ring_spectrum_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> F, 
	Array<DP, 1> sector_angle_array, 
	DP n, 
	Array<DP,2> local_Sk, 
	DP kfactor[]
);
								
/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.								
void Compute_ring_spectrum_SCFT
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
void Compute_local_ring_spectrum_SCFT
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
void Compute_ring_spectrum_SCFT
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


																				
/** @brief Computes local shell spectrum \f$ \sum \vec{K'} . (\vec{Vr}(K') 
*					\times \vec{Vi}(K'))) \f$ in local processor.
 *
 * @param  N[]  The size of the array A.
 * @param  Ax  x-component of vector field \f$ \vec{A} \f$.
 * @param  Ay  y-component of the vector field \f$ \vec{A} \f$.
 * @param  Az  z-component of the vector field \f$ \vec{A} \f$.
 * @param  kfactor[]  factor to compute actual wavenumber given grid wavenumber.
 *
 * @return H1(K) = \f$ \sum \vec{K'} . (\vec{Vr}(K') \times \vec{Vi}(K'))) 
 *						\f$ for  \f$ K \le K' K+1 \f$ (K>0).
 */												 
void Compute_local_ring_spectrum_helicity_SCFT
(
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	Array<DP, 1> sector_angle_array,
	Array<DP,2> local_H1k, 
	DP kfactor[]
);	

/// Returns sum all local_E's.
/// Remove the origin's contribution if master processor.													
void Compute_ring_spectrum_helicity_SCFT
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
void Compute_local_cylinder_ring_spectrum_SCFT
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
void Compute_cylinder_ring_spectrum_SCFT
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
void Compute_local_cylinder_ring_spectrum_SCFT
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
void Compute_cylinder_ring_spectrum_SCFT
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
void Compute_local_cylinder_ring_spectrum_SCFT
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
void Compute_cylinder_ring_spectrum_SCFT
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
void Compute_local_cylinder_ring_spectrum_SCFT
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
void Compute_cylinder_ring_spectrum_SCFT
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
void Compute_local_cylinder_ring_spectrum_helicity_SCFT
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
void Compute_cylinder_ring_spectrum_helicity_SCFT
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
void Compute_local_imag_shell_spectrum_B0_SCFT
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
void Compute_imag_shell_spectrum_B0_SCFT
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
void Compute_local_imag_ring_spectrum_B0_SCFT
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
void Compute_imag_ring_spectrum_B0_SCFT
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
void Compute_local_imag_cylinder_ring_spectrum_B0_SCFT
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
void Compute_imag_cylinder_ring_spectrum_B0_SCFT
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

//*****************************  End of scft_energy.h *****************************************	














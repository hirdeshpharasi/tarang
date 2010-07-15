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

/*! \file  Csf.h
 * 
 * @brief  Class declaration of Csf, a Complex Scalar Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bugs Out of date functions CS_output_split_arrays, CS_iutput_split_arrays.
 */

#ifndef _CSF
#define _CSF

#include "array_fn/array_basic_inline.h"
#include "array_fn/array_basic.h"
#include "array_fn/struct_fn.h"
#include "array_fn/planar_struct_fn.h"

#include "universal/universal_inline.h" 
#include "universal/universal_basic.h" 
#include "universal/universal_tr.h" 
#include "universal/universal_energy.h" 
#include "universal/universal_ET.h" 
#include "universal/universal_misc.h" 

//*********************************************************************************************

//! Complex Scalar field
/*!
 * 3 dimensional complex scalar field with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, N_3/2+1] \f$.   <BR>
 * This array can also store real values, and its dimension is \f$[N_1, N_2, N_3] \f$. 
 *
 * Inverse transform: Complex -> Real.  <BR>
 * Forward transform: Real -> Complex. <BR>
 * The implemented basis functions: <BR>
 *  FOUR: Fourier along all directions. <BR>
 *  SCFT: Sin/Cos along x, and Fourier along all perpendicular directions. <BR>
 * 
 *  The modal energy \f$ C(\vec{k}) \f$ is defined in each basis function.  
 *  The dissipation rate for the mode is \f$ 2 K^2 C(\vec{k}) \f$.
 * 
 *  Entropy = = \f$ \sum p(k) log(1/p(k)) \f$ where probability  \f$ p(k) =  E(k)/E_{total)\f$ 
 *					with \f$ E(k) = |A(k)|^2/2 \f$.
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */



class CSF 
{ 

public:

	// DYNAMIC ARRAYS
	
	//!  \f$F(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *F;				
	
	//!  Energy of F in shell k								
	Array<DP,1>		*CS_shell_ek;					
	
	//!  Energy Dissipation rate in shell k (without \f$ \kappa \f$).
	Array<DP,1>		*CS_shell_dissk;				
	
	//!  Energy in ring(m,n)
	Array<DP,2>		*CS_ring_ek;					
	
	//!  Energy dissipation in ring(m,n) (without \f$ \kappa \f$).
	Array<DP,2>		*CS_ring_dissk;		
	
	//!  Energy in ring(m,n) 
	Array<DP,2>		*CS_cylinder_ring_ek;					
	
	//!  Energy Dissipation rate in ring(m,n) (without \f$ \nu \f$).
	Array<DP,2>		*CS_cylinder_ring_dissk;			


	// OTHER VARS
	//!  Basis type: FOUR or SCFT
	string	CS_basis_type;	
	
	//!  Alias_switch: ALIAS or DEALIAS.
	string  CS_alias_switch;
	
	//! Number input field mode: ASCII, BINARY
	string CS_no_input_field_mode;
	
	//! Number input field mode: ASCII, BINARY
	string CS_no_output_field_mode;
	
	//! 1 if anisotropy_ring is on; 0 otherwise
	int		CS_anisotropic_ring_switch;
	
	//! 1 if anisotropy_cylinder is on; 0 otherwise
	int		CS_anisotropic_cylinder_switch;
	
	//! 1 if on; 0 if off
	int		CS_structure_fn_switch;
	
	//! 1 if on; 0 if off
	int		CS_planar_structure_fn_switch;
	
	//!  Size of array
	int		Ncs[4];	
	
	//!  Conversion factor for grid to actual wavenumber: \f$ f_i \f$.																										
	DP		CS_kfactor[4];	
	
	//!  Conversion factor for grid to actual real space position: \f$ f_i \f$.
	DP		CS_xfactor[4];
		
	//!  Total energy of Csf.													
	DP		CS_total_energy;						
	
	//!  Total dissipation rate of Csf (without \f$ \kappa \f$).
	DP		CS_total_dissipation;					
	
	//!  Entropy of Csf
	DP		CS_entropy;																				
	
		
	
	//!  Size of CS_shell_ek	
	int		CS_shell_ek_size;	
	
	//!  Size of CS_ring_ek.	
	int		CS_ring_ek_size;
	
	//!  Scheme = 0: equispaced angles; scheme=1; angles s.t. each sector has equal modes.
	int		CS_sector_spectrum_input_scheme;										
	
	//!  Number of sectors for ring spectrum.
	int		CS_no_sectors_spectrum;					
	
	//!	ring_n = [ sector_angle(n-1),sector_angle(n) );
	/// the last ring includes theta=pi.
	Array<DP,1>  *CS_sector_angle_array_spectrum;			
													

	//! Size of CV_cylinder_shell_ek_size
	int	CS_cylinder_shell_ek_size;
	
	//! Number of slabs in cylinder for spectrum
	int CS_no_cylinder_slabs_spectrum;
	
	//! Cylinder Kpll array for spectrum; ring(m) = [H(m-1), H(m)].
	Array<DP,1> *CS_cylinder_kpll_array_spectrum;
	
	
	// min q for structure function
	int	CS_structurefn_q_min;
	
	// max q for structure function
	int	CS_structurefn_q_max;
	
	//! Max r for structure function
	int	CS_structurefn_r_max;
	
	//! max grid coordinate along the anisotropy direction.
	int CS_structure_fn_rpll_max;
	
	//! structure function
	/*
	 *	St(r, q, m) with  r=[0, rmax], q=[q_min, q_max], m=0:1
	 *  rmax = maximum radius that fits inside real field.
	 *  m=0:  for \f$ \Delta u_{||} \f$
	 *  m=1:  for \f$ \Delta u_{\perp} \f$
	 */
	Array<DP,2> *CS_St;
	
	Array<DP,3> *CS_st_planar;
												
												
//*********************************************************************************************					
	
public:

	/*! A constructor; 
	 *
	 *  Allocation for F, Ek etc. 
	 * Initialize the arrays.
	 */
	CSF
	(
		int *NN, 
		string string_switches[], 
		Array<int,1> switches,
		DP *kfactor, 
		Array<int,1> misc_output_para
	);


	/// Set F=0.
	void CS_Initialize();
	
	
	/*! @brief Inplace Forward transform of CSF; \f$ \mathcal{F}(F) -> F \f$.
	*
	*  @param	 temp_r			Complex 3D array, a temporary array.
	*
	*  @return FOUR: FOURIER transform
	*  @return SCFT: SFT(F) -> F.
	*/
	void CS_Forward_transform(Array<complx,3> temp_r);
	
	
	/*! @brief Inplace Inverse transform of CSF; \f$ \mathcal{F}^{-1}(F) -> F \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(F) -> F.
	 */	
	void CS_Inverse_transform(Array<complx,3> temp_r);
	
	
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}(F) -> F	\f$.
	 *
	 *  @param	 Ftr			Ftr that is to be transformed. 
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(Ftr) -> F  (Ftr unaffected).
	 */	
	void CS_Forward_transform_transpose_order(Array<complx,3> Ftr, Array<complx,3> temp_r);
	
	
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}^{-1}(F) -> Ftr	\f$.
	 *
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(F) -> Ftr  (F unaffected).
	 */		
	void CS_Inverse_transform_transpose_order(Array<complx,3> Ftr, Array<complx,3> temp);


	//***************************************************************************************** 
	
	/// F(k) = F(k)/K^2.  F(0)=0.
	void CS_divide_ksqr();										
  	
	//***************************************************************************************** 
	
	/// Compute total energy and dissipation of F.
	void CS_Compute_totalenergy_diss();
	
	/// Compute entropy of F.
	void CS_Compute_entropy();
	
	//***************************************************************************************** 
	
	/// Compute shell spectrum of F.
	void CS_Compute_shell_spectrum();
	
	/// Compute shell spectrum C(F,G).							
	void CS_Compute_shell_spectrum(Array<complx,3> G, Array<DP,1> result);
	
	
	/// Compute ring spectrum of \f$ \vec{V} \f$.
	void CS_Compute_ring_spectrum();
	
	void CS_Compute_ring_spectrum(Array<complx,3> G, Array<DP,2> result);
	
	/// Compute cylinderical ring spectrum of F.
	void CS_Compute_cylinder_ring_spectrum();
	
	void CS_Compute_cylinder_ring_spectrum(Array<complx,3> G, Array<DP,2> result);
	
	//*****************************************************************************************	
	
	/// 3D: Return modal energy for grid index (i1,i2,i3).
	DP CS_Modal_energy(int i1, int i2, int i3);
	
	//*****************************************************************************************	
	
	/// Dealiase F.
	void CS_Dealias();
	
	//*****************************************************************************************	
	
	/// Output Csf F to file_out.	temp_array is a temporary array.			
	void CS_output(ofstream& file_out, Array<complx,3> temp_array);
	
	/// Output reduced F to file_out.	temp_array is a temporary array.
	void CS_output(ofstream& file_out, int Nreduced[], Array<complx,3> temp_array);
	
	/// Input F from file_in.  temp_array is a temporary array.
	void CS_input(ifstream& file_in, Array<complx,3> temp_array);
	
	/// Input reduced F from file_in.  temp_array is a temporary array.
	void CS_input(ifstream& file_in, int Nreduced[], Array<complx,3> temp_array);
	
	
}; 

#endif

//******************************** End of CSF.h ***********************************************




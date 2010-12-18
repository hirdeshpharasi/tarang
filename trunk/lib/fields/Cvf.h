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

/*! \file  Cvf.h
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008 
 * @bugs  No known bugs
 */

#ifndef _CVF
#define _CVF

#include "universal/universal_inline.h" 
#include "universal/universal_basic.h" 
#include "universal/universal_tr.h" 
#include "universal/universal_energy.h" 
#include "universal/universal_ET.h"  
#include "universal/universal_misc.h" 

#include "array_fn/array_basic_inline.h"
#include "array_fn/array_basic.h"
#include "array_fn/struct_fn.h"
#include "array_fn/planar_struct_fn.h"


//*********************************************************************************************

//! Complex vector field
/*!
 * 3 dimensional complex vector field has 3 components with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2,N_3/2+1] \f$.   <BR>
 * These arrays can also store real values, and their dimension are \f$[N_1, N_2, N_3] \f$. 
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
 *  Helicity1 = \f$ H1(K) = \vec{K} . [\vec{Vr} \times \vec{Vi}] \f$. <BR>
 *  Helicity2 = \f$ H2(K) = \vec{K} . [(\vec{Vr} \times \vec{Vi})] /K^2 \f$. <BR>
 * 
 *  Entropy = \f$ \sum p(k) log(1/p(k)) \f$ where probability  
 *				\f$ p(k) =  E(k)/E_{total) \f$ with \f$ E(k) = |A(k)|^2/2 \f$.
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */

//*********************************************************************************************

class CVF 
{
 
public:
	//!  \f$ V_x(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *V1;							
	
	//!  \f$ V_y(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *V2;							
	
	//!  \f$ V_z(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *V3;								

	//!  Energy of Vx in shell k
	Array<DP,1>		*CV_shell_ek1;					
	
	//!  Energy of Vy in shell k
	Array<DP,1>		*CV_shell_ek2;					
	
	//!  Energy of Vz in shell k
	Array<DP,1>		*CV_shell_ek3;
									
	//!  Sum of the above energies of Vx, Vy, and Vz in shell k
	Array<DP,1>		*CV_shell_ek;					
	
	//!  Energy Dissipation rate in shell k (without \f$ \nu \f$).
	Array<DP,1>		*CV_shell_dissk;										
	
	//! Sum \f$ \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K})) \f$ in shell k
	//! components along 1,2,3 directions.
	Array<DP,1>		*CV_shell_h1k1;					
	Array<DP,1>		*CV_shell_h1k2;	
	Array<DP,1>		*CV_shell_h1k3;							
	
	//!  Energy in ring(m,n) along e1 (toroidal direction).
	Array<DP,2>		*CV_ring_ek1;					
	
	//!  Energy in ring(m,n) along e2 (poloidal direction).
	Array<DP,2>		*CV_ring_ek2;					
	
	//!  Energy Dissipation rate in ring(m,n) (without \f$ \nu \f$).
	Array<DP,2>		*CV_ring_dissk;						
	
	//!   \f$ \sum \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K})) \f$ 
	/// in ring(m,n)
	Array<DP,2>		*CV_ring_h1k;					


	//!  Energy spectrum along the anisotropy direction.
	Array<DP,2>		*CV_cylinder_ring_ek1;					
	
	//!  Energy spectrum perpendicular to the anisotropy direction.
	Array<DP,2>		*CV_cylinder_ring_ek2;					
	
	//!  Energy Dissipation rate in ring(m,n) (without \f$ \nu \f$).
	Array<DP,2>		*CV_cylinder_ring_dissk;					
	
	//! \f$ \sum \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K})) \f$ 
	/// in ring(m,n)
	Array<DP,2>		*CV_cylinder_ring_h1k;	
	
	//*****************************************************************************************
	
	// Nondynamic allocation now
	
	//!  Basis type: FOUR or SCFT
	string  CV_basis_type;	
	
	//!  Alias_switch: ALIAS or DEALIAS.
	string  CV_alias_switch;
	
	string CV_no_input_field_mode;
	string CV_no_output_field_mode;
	
	//! 1 if anisotropy_ring is on; 0 otherwise
	int		CV_anisotropic_ring_switch;
	
	//! 1 if anisotropy_cylinder is on; 0 otherwise
	int		CV_anisotropic_cylinder_switch;
	
	//! 1 if on; 0 if off
	int		CV_structure_fn_switch;
	
	//! 1 if on; 0 if off
	int		CV_planar_structure_fn_switch;
	
	//!  Size of array							
	int		Ncv[4];	
	
	//!  Conversion factor for grid to actual wavenumber: \f$ f_i \f$.																		
	DP		CV_kfactor[4];							
	
	//!  Conversion factor for grid to actual real space position: \f$ f_i \f$.
	DP		CV_xfactor[4];
	
	//!  Total energy of Cvf.
	DP		CV_total_energy;						
	
	//!  Total dissipation rate of Cvf (without \f$ \nu \f$).
	DP		CV_total_dissipation;	
	
	//!  \f$ \sum \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K})) \f$				
	DP		CV_total_helicity1;						
	
	//!  \f$ \sum \vec{K} \cdot (\Re\vec{V}(\vec{K}) \times  \Im\vec{V}(\vec{K}))/K^2 \f$
	DP		CV_total_helicity2;						
	
	//!  Dissipation of H1. 
	DP		CV_total_dissipation_H1;				
	
	//!  Dissipation of H2. 
	DP		CV_total_dissipation_H2;				
	
	//!  Entropy of Cvf
	DP		CV_entropy;	
	
	//!  Size of CV_shell_ek.
	int		CV_shell_ek_size;
	
	//!  Size of CV_ring_ek
	int		CV_ring_ek_size;						
	
	//!  Scheme = 0: equispaced angles; scheme=1; angles s.t. each sector has equal modes.
	int		CV_sector_spectrum_input_scheme;
	
	//!  Number of sectors for ring spectrum.
	int		CV_no_sectors_spectrum;	
					
	//!  sector_angles array; ring_n = ( sector_angle(n-1),sector_angle(n) ].
	// the first ring includes theta=0.
	Array<DP,1>  *CV_sector_angle_array_spectrum;			
													
		
	//! Size of CV_cylinder_shell_ek_size
	int	CV_cylinder_shell_ek_size;
	
	//! Number of slabs in cylinder for spectrum
	int CV_no_cylinder_slabs_spectrum;
	
	//! Cylinder Kpll array for spectrum; ring(m) = [H(m-1), H(m)].
	Array<DP,1> *CV_cylinder_kpll_array_spectrum;
	
	
	// min q for structure function
	int	CV_structurefn_q_min;
	
	// max q for structure function
	int	CV_structurefn_q_max;
	
	//! Max r for structure function
	int	CV_structurefn_r_max;
	
	//! max grid coordinate along the anisotropy direction.
	int CV_structure_fn_rpll_max;
	
	//! structure function
	/*
	 *	CV_St(r, q, m) with  r=[0, rmax], q=[q_min, q_max], m=0:1
	 *  rmax = maximum radius that fits inside real field.
	 *  m=0:  for \f$ \Delta u_{||} \f$
	 *  m=1:  for \f$ \Delta u_{\perp} \f$
	 */
	Array<DP,3> *CV_St;
	
	//! Planar structure function
	/*
	 *	CV_St_planar(kpll, r, q, m) with  r=[0, rmax], q=[q_min, q_max], m=0:1
	 *  rmax = maximum radius that fits inside real field.
	 *  m=0:  for \f$ \Delta u_{||} \f$
	 *  m=1:  for \f$ \Delta u_{\perp} \f$
	 */
	Array<DP,4> *CV_st_planar;				
																						
	//*****************************************************************************************
	
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */
	CVF
	(
		int *NN, 
		string string_switches[], 
		Array<int,1> switches,
		DP *kfactor, 
		Array<int,1> misc_output_para
	);
 
	/// Copy curr Vi's to VF to.
	void CV_Copy_to(CVF& to);  
	
	/// Set \f$ \vec{V} = 0 \f$.  
	void CV_Initialize();
	
	// Initialize fftw_plans
	void Init_fftw_plan();
	

	//*****************************************************************************************	
	
	/*! @brief Inplace Forward transform of CVF; \f$ \mathcal{F}(\vec{V}) -> \vec{V} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: FOURIER transform
	 *  @return SCFT: SFT(V1) -> V1; CFT(V2) -> V2; CFT(V3) -> V3;
	 */
	void CV_Forward_transform(Array<complx,3> temp_r);
	

	/*! @brief  Forward transform in transpose order of CVF; 
	 *			\f$ \mathcal{F}(\vec{Vtr}) -> \vec{V}	\f$.
	 *
	 *  @param	 V1tr			x component of Vtr that is to be transformed. 
	 *  @param	 V2tr			y component of Vtr that is to be transformed. 
	 *  @param	 V3tr			z component of Vtr that is to be transformed. 
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(V1tr) -> V1; CFT(V2tr) -> V2; CFT(V3tr) -> V3;  (Vtr unaffected).
	 */	
	void CV_Forward_transform_transpose_order
	(
		Array<complx,3> V1tr, 
		Array<complx,3> V2tr,
		Array<complx,3> V3tr, 
		Array<complx,3> temp_r
	);
													
	//*****************************************************************************************													
													
	/*! @brief Inplace Inverse transform of CVF; \f$ \mathcal{F}(\vec{V}) -> \vec{V} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(V1) -> V1; ICFT(V2) -> V2; ICFT(V3) -> V3;
	 */													
	void CV_Inverse_transform(Array<complx,3> temp_r);		
	
	
	
	/*! @brief  Inverse transform in transpose order of CVF; 
	 *			\f$ \mathcal{F}(\vec{V}) -> \vec{Vtr}	\f$.
	 *
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(V1) -> V1tr.
	 *  @return SCFT: ICFT(V2) -> V2tr.
	 *  @return SCFT: ICFT(V3) -> V3tr.
	 */																								
	void CV_Inverse_transform_transpose_order
	(
		Array<complx,3> V1tr, 
		Array<complx,3> V2tr,
		Array<complx,3> V3tr, 
		Array<complx,3> temp
	);	
													
	 
	//*****************************************************************************************   
	
	/// Compute total energy and dissipation of \f$ \vec{V} \f$. 
	void CV_Compute_totalenergy_diss();
	
	/// Compute entropy of \f$ \vec{V} \f$.
	void CV_Compute_entropy();
	
	//*****************************************************************************************	
	
	/// Compute shell spectrum of \f$ \vec{V} \f$.
	void CV_Compute_shell_spectrum();

	/// Compute shell spectrum \f$ C(\vec{V}, \vec{W}) \f$.
	void CV_Compute_shell_spectrum
	(
		Array<complx,3> W1, Array<complx,3> W2, Array<complx,3> W3, 
		Array<DP,1> result
	);	
	
	/// Compute ring spectrum of \f$ \vec{V} \f$.
	void CV_Compute_ring_spectrum();
	
	/// Compute ring spectrum \f$ C(\vec{V}, \vec{W}) \f$.
	void CV_Compute_ring_spectrum
	(
		Array<complx,3> W1, Array<complx,3> W2, Array<complx,3> W3,
		Array<DP,2> result
	);
	
	// Compute cylinderical ring spectrum of \f$ \vec{V} \f$.
	void CV_Compute_cylinder_ring_spectrum();
	
	/// 3D: Compute cylinderical ring spectrum \f$ C(\vec{V}, \vec{W}) \f$.
	void CV_Compute_cylinder_ring_spectrum
	(
		Array<complx,3> W1, Array<complx,3> W2, Array<complx,3> W3, 
		Array<DP,2> result
	);
	
	//*****************************************************************************************
	
	
	/// Return modal energy for grid index (i1,i2,i3).
	DP CV_Modal_energy(int i1, int i2, int i3);
	
	/// vorticity = modal vorticity for grid index (i1,i2,i3).
	void CV_Modal_vorticity(int i1, int i2, int i3, TinyVector<complx,3> &vorticity);
	
	//*****************************************************************************************
	
	/// Return modal helicity for the grid index (i1,i2,i3).
	DP CV_Modal_helicity(int i1, int i2, int i3);
	
	/// Compute total helicity1, helicity2, and their dissipation for Cvf.
	void CV_Compute_total_helicity();
	
	/// Compute helicity spectrum (H1) of Cvf for the shells.
	void CV_Compute_shell_spectrum_helicity();	
	
	/// Compute helicity spectrum (H1) of Cvf for the rings.
	void CV_Compute_ring_spectrum_helicity();
	
	/// Compute helicity spectrum (H1) of Cvf for the cylindrical rings.
	void CV_Compute_cylinder_ring_spectrum_helicity();			
	
	//*****************************************************************************************	
	
	// Compute structure function
	void CV_Compute_structure_function();
	
	void CV_Compute_planar_structure_function();
	
	
	//*****************************************************************************************	
	
	 /// Dealiase \f$ \vec{V} \f$ using 2/3 rule.
	 void CV_Dealias();

	//*****************************************************************************************																																			
	 
	/// Output \f$ \vec{V} \f$  to file_out.	temp_array is a temporary array.																																			
	void CV_output(ofstream& file_out, Array<complx,3> temp_array, string format);
	
//	void CV_output_hdf5(DataSet* dataset1R, DataSet* dataset1C, DataSet* dataset2R,
//			DataSet* dataset2C, DataSet* dataset3R, DataSet* dataset3C,
//		       	Array<complx,3> temp_array);

	/// Output reduced \f$ \vec{V} \f$  to file_out.	
	void CV_output(ofstream& file_out, int Nreduced[], string format);
	
	/// Input \f$ \vec{V} \f$ from file_in.  temp_array is a temporary array.
	void CV_input(ifstream& file_in, Array<complx,3> temp_array, string format);
	
	/// Input reduced \f$ \vec{V} \f$ from file_in.  
	void CV_input(ifstream& file_in, int Nreduced[], string format);
	
	//*****************************************************************************************
	
};

#endif

//******************************** End of CVF.h ***********************************************	







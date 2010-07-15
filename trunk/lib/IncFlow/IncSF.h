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

/*! \file  IncSF.h
 * 
 * @brief  Class declaration of IncSF, Incompressible scalar Field 
 *		example: passive scalar, temperature in RB convection
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug  No known bugs
 */

//*********************************************************************************************

#ifndef _IncSF
#define _IncSF

#include "../fields/fields.h"
#include "../plainfields/plainfields.h"   
#include "Nlin.h"

//! @brief Incompressible scalar field IncSF 
/*!
 *  Inherits CSF that contains the complex scalar field. <BR>
 *  RSF that contains real vector field, typically Inverse tranform of CSF. <BR>
 * 
 *	Compute nonlinear terms <BR>
 *  Compute energy transfer functions: <BR>
 *	Isotropic: flux, shell-to-shell <BR>
 *  Anisotropic: ring-to-ring in spherical shell and in cylinderical shells.
 *
 *	@sa IncSF.h
 *	@sa Nlin.h
 *  @sa EnergyTr.h
 */
 
 
//*********************************************************************************************

class IncSF: public CSF, public RSF 
{ 
  public:
  
	//! nlin = \f$ N^\psi = \vec{V} \cdot \nabla \psi \f$.
	Array<complx,3> *nlin;  
	
	//!  temp array  \f$ (local_{N1}, N_2, N_3/2+1) \f$.       
	Array<complx,3> *SF_temp;
	
	//!  Force array for scalar \f$ (local_{N1}, N_2, N_3/2+1) \f$.
	Array<complx,3> *Force;  
	
	//!  spectrum of energy feed by force 
	/// \f$ \sum \Re(F(\vec{K}) \cdot \psi^*(\vec{K})) \f$.
	Array<DP,1>		*shell_spectrum_force_SFk;		 

	//! Size of complex array of CSF, force etc.
	int		NIs[4];		
					
	string	ISF_basis_type;
	string	ISF_alias_switch;
	
	DP		ISF_kfactor[4];
	
	//! diffusivity of the scalar
	DP		diffusion_coefficient;			
	int		hyper_diffusion_switch;
	DP		hyper_diffusion_coefficient;
	
	// Anisotropic arrays
	Array<DP,2>		*ring_spectrum_force_SFk;
	Array<DP,2>		*cylinder_ring_spectrum_force_SFk;
	
//*********************************************************************************************	
	
 public:
	
	IncSF
	(
		int NN[], 
		string string_switches[], 
		Array<int,1> switches, 
		DP *prog_kfactor, 
		DP diff_coefficient,  
		DP hyper_diff_coefficient, 
		Array<int,1> misc_output_para
	);
	
	void Mult_field_exp_ksqr_dt(DP dt);
	void Mult_nlin_exp_ksqr_dt(DP dt);
	
	void Copy_field_to(CSF& T);
	void Copy_field_to(PlainCSF& T);
	void Copy_field_from(CSF& T);
	void Copy_field_from(PlainCSF& T);
	
	void Add_nlin_dt(DP dt);
	void Add_nlin_to_field(CSF& T, DP factor);
	void Add_nlin_to_field(PlainCSF& T, DP factor);
	
	void Compute_force_shell_spectrum();
	void Compute_force_ring_spectrum();
	void Compute_force_cylinder_ring_spectrum();	
	
	void Add_complex_conj(int kx, int ky, int kz, complx G);
	void Add_complex_conj_force(int kx, int ky, int kz, complx ForceG);
	
	void Assign_field_comp_conj(int kx, int ky, int kz, complx G);
	void Assign_force_and_comp_conj(int kx, int ky, int kz, complx FG);
	
	DP Get_Tk(int kx, int ky, int kz);
	
	void Put_scalar_amp_phase_comp_conj(int i1, int i2, int i3, DP amp, DP phase);
	
	void Force_scalar_alpha(int kx, int ky, int kz, DP alpha);
};

#endif

//************************** End of IncSF.h   *************************************************	




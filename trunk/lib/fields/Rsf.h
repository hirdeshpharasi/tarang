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

/*! \file  Rsf.h
 * 
 * @brief  Class declaration of Rsf, a Real Scalar Field 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * 
 * @bug  No known bugs
 */
 

#ifndef _RSF
#define _RSF

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

//! Real Scalar field
/*!
 * 3 dimensional real scalar field with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, N_3/2+1] \f$.   <BR>
 * This array contains two real values for every complex value.
 *  The dimension of the real array is \f$[N_1, N_2, ...N_3] \f$. 
 *
 * Inverse transform: Complex -> Real.  <BR>
 * Forward transform: Real -> Complex. <BR>
 * The implemented basis functions: <BR>
 *  FOUR: Fourier along all directions. <BR>
 *  SCFT: Sin/Cos along x, and Fourier along all perpendicular directions. <BR>
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */

//*********************************************************************************************

class RSF 
{ 
public:

	//!  \f$F(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *Fr;
	
	//!  Size of array
	int		Nrs[4];								
	
	//!  Basis type: FOUR or SCFT
	string	RS_basis_type;	

//*********************************************************************************************
public:

	/*! A constructor; 
	 *
	 *  Allocation for Fr, Ek etc. 
	 * Initialize the arrays.
	 */
	RSF(int NN[], string string_switches[]);

	//*****************************************************************************************
	
	/*! @brief Inplace Forward transform of CSF; \f$ \mathcal{F}(Fr) -> Fr \f$.
	*
	*  @param	 temp_r			Complex 3D array, a temporary array.
	*
	*  @return FOUR: FOURIER transform
	*  @return SCFT: SFT(Fr) -> Fr.
	*/
	void RS_Forward_transform(Array<complx,3> temp_r);
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}(Fr) -> F	\f$.
	 *
	 *  @param	 Fr				Fr that is to be transformed. 
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(Fr) -> F  (Ftr unaffected).
	 */	
	void RS_Forward_transform_transpose_order(Array<complx,3> F, Array<complx,3> temp_r);
	
	
	//*****************************************************************************************
	
	/*! @brief Inplace Inverse transform of CSF; \f$ \mathcal{F}^{-1}(Fr) -> Fr \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(Fr) -> Fr.
	 */	
	void RS_Inverse_transform(Array<complx,3> temp_r);
	
	
	/*! @brief  Forward transform in transpose order of CSF; \f$ \mathcal{F}^{-1}(F) -> Fr	\f$.
	 *
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(F) -> Fr  (F unaffected).
	 */	
	void RS_Inverse_transform_transpose_order(Array<complx,3> F, Array<complx,3> temp);

	//*****************************************************************************************
	
	/// dealias F.
	void RS_Dealias();
	
	//*****************************************************************************************
	
	/// Output F to fileout as real.
	void RS_Output(ofstream& fileout, Array<complx,3> temp_array);
	
}; 

#endif

//**************************  End of RSF class declaration ************************************




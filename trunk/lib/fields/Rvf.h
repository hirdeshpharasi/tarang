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

/*! \file  Rvf.h
 * 
 * @brief  Class declaration of Rvf, a Real Vector Field 
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug  No known bugs
 */
 
 
#ifndef _RVF
#define _RVF

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

/*

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

 */ 
//*********************************************************************************************

//! Real vector field
/*!
 * 3 dimensional real vector field has 3 components with 3 array indices. <BR>
 * The size of the arrays are \f$[N_1, N_2, N_3/2+1] \f$.   <BR>
 * These arrays hold two real values for each complex entry.
 *  The real dimension are \f$[N_1, N_2, N_3] \f$. 
 *
 * Inverse transform: Complex -> Real.  <BR>
 * Forward transform: Real -> Complex. <BR>
 * The implemented basis functions: <BR>
 *  FOUR: Fourier along all directions. <BR>
 *  SCFT: Sin/Cos along x, and Fourier along all perpendicular directions. <BR>
 * 
 *
 *	@sa field_intermal_four.h
 *	@sa field_internal_sincosfour.h
 */

//*********************************************************************************************

class RVF 
{ 

public:
	//!  \f$ V_x(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *V1r;	
	
	//!  \f$ V_y(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *V2r;			
	
	//!  \f$ V_z(local_N1, N_2, N_3/2+1) \f$.
	Array<complx,3> *V3r;		

	//!  Size of array
	int		Nrv[4];	
	
	//!  Basis type: FOUR or SCFT						
	string	RV_basis_type;	
	
//*********************************************************************************************			
 
public:

	/*! A constructor; 
	 *
	 *  Allocation for Vi, Ek etc. 
	 * Initialize the arrays.
	 */
	RVF(int *NN, string string_switches[]);

	//*****************************************************************************************
	
	/*! @brief Inplace Forward transform of CVF; \f$ \mathcal{F}(\vec{Vr}) -> \vec{Vr} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: FOURIER transform
	 *  @return SCFT: SFT(V1r) -> V1r; CFT(V2r) -> V2r; CFT(V3r) -> V3r;
	 */
	void RV_Forward_transform(Array<complx,3> temp_r);
	
	//
	//
	
	/*! @brief  Forward transform in transpose order of CVF; 
	 *				\f$ \mathcal{F}(\vec{Vr}) -> \vec{V}	\f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: SFT(V1r) -> V1; CFT(V2r) -> V2; CFT(V3r) -> V3;  (Vr unaffected).
	 */	
	void RV_Forward_transform_transpose_order
	(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp_r
	);
		
	//*****************************************************************************************
																									
	/*! @brief Inplace Inverse transform of CVF; \f$ \mathcal{F}(\vec{Vr}) -> \vec{Vr} \f$.
	 *
	 *  @param	 temp_r			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Inverse FOURIER transform
	 *  @return SCFT: ISFT(V1r) -> V1r; ICFT(V2r) -> V2r; ICFT(V3r) -> V3r;
	 */													
	void RV_Inverse_transform(Array<complx,3> temp_r);

	//
	//
													
	/*! @brief  Inverse transform in transpose order of CVF; 
	 *				\f$ \mathcal{F}(\vec{V}) -> \vec{Vr}	\f$.
	 *
	 *  @param	 V1				x component of V that is to be transformed. 
	 *  @param	 V2				y component of V that is to be transformed. 
	 *  @param	 V3r			z component of V that is to be transformed. 
	 *  @param	 temp			Complex 3D array, a temporary array.
	 *
	 *  @return FOUR: Not applicable so far.
	 *  @return SCFT: ISFT(V1) -> V1r.
	 *  @return SCFT: ICFT(V2) -> V2r.
	 *  @return SCFT: ICFT(V3) -> V3r.
	 */																
	void RV_Inverse_transform_transpose_order
	(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp
	);
												
	//*****************************************************************************************
												
	/*! @brief Inplace Forward transform of real space product Vi*Vj.
	  *
	  *	Before the Forward transform, the state Vr is <BR>
	  * V1r = V1r*V2r (Odd parity); V2r = V2r*V3r (Even parity); V3r = V3r*V1r (Odd parity)	
	  *  
	  *
	  *  @param	 temp_r			Complex 3D array, a temporary array (transposed order).
	  *
	  *  @return  \f$ V1r = \mathcal{F}(V1r*V2r);  V2r = \mathcal{F}(V2r*V3r); 
	  *				V3r = \mathcal{F}(V3r*V1r) \f$.
	  */																																		
	void RV_Forward_transform_RSprod(Array<complx,3> temp_r);
	
	
	
	/*! @brief Transpose_order Forward transform of real space product Vi*Vj.
	  *
	  *	Before the Forward transform, the state Vr is <BR>
	  * V1r = V1r*V2r (Odd parity); V2r = V2r*V3r (Even parity); V3r = V3r*V1r (Odd parity). 	 
	  *
	  *  @param	 temp_r			Complex 3D array, a temporary array (transposed order).
	  *
	  *  @return  \f$ V1 = \mathcal{F}(V1r*V2r);  V2 = \mathcal{F}(V2r*V3r); 
	  *					V3 = \mathcal{F}(V3r*V1r) \f$.
	  */
	void RV_Forward_transform_RSprod_transpose_order
	(
		Array<complx,3> V1, Array<complx,3> V2, Array<complx,3> V3, 
		Array<complx,3> temp_r
	);																								
	
	
	//*****************************************************************************************
	/// Output \f$ \vec{Vr} \f$ as real to fileout.
	void RV_Output(ofstream& fileout, Array<complx,3> temp_array, string format);
	
	void RV_Output_transpose_order(ofstream& fileout, Array<complx,3> temp_array, string format);
	
	//void RV_Output_hdf5(DataSet* dataset1, DataSet* dataset2, DataSet* dataset3,
	//	       	 Array<complx,3> temp_array);
	
	//*****************************************************************************************
	/// Dealias \f$ \vec{Vr} \f$  using 2/3 rule.
	void RV_Dealias();
};

#endif

//**************************  End of RVF.h  ***************************************************



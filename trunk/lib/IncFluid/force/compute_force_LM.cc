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

/*! \file  compute_force_LM.cc
 * 
 * @brief  Set up force using Liquid metal for under strong mean magnetic field.
 *
 * @note F =  (hat(B0).hat(K))^2 V
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"


//*********************************************************************************************


void IncFluid::Compute_force_Liquid_metal()
{	
	
	DP B0x = (*force_double_para)(1);
	DP B0y = (*force_double_para)(2);
	DP B0z = (*force_double_para)(3);
	
	TinyVector <DP,3> B0;
	B0 = B0x, B0y, B0z;
	
	*Force1 = *V1;	
	*Force2 = *V2; 
	*Force3 = *V3;
	
	Array_mult_V0_khat_sqr(basis_type, N, *Force1, B0, kfactor);
	Array_mult_V0_khat_sqr(basis_type, N, *Force2, B0, kfactor);
	Array_mult_V0_khat_sqr(basis_type, N, *Force3, B0, kfactor);
	
	*Force1 = complx(-1,0) * (*Force1);
	*Force2 = complx(-1,0) * (*Force2);
	*Force3 = complx(-1,0) * (*Force3);
	
	if (alias_switch == "DEALIAS")   Dealias_force();

}	

//****************************** End of compute_force_LM.cc ***********************************



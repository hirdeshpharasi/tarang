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


/*! \file  compute_force_TG.cc
 * 
 * @brief  Set up force using Taylor Green theory.
 *
 * @note Fx = F*sin(k0 x) cos(k0 y) cos(k0 z) <BR>
 *		Fy = -F*cos(k0 x)sin(k0 y) cos(k0 z) <BR>
 *		Fz = 0
 *
 *  @note is_force_field_para_read is used to read only once.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"

//*********************************************************************************************


void IncFluid::Compute_force_Taylor_Green()
{

	if (is_force_field_para_read == 0)
	{
		int k0 = ((int) (*force_field_para)(1));
		DP force_amp = (*force_field_para)(2);
		
		Setup_Taylor_Green_force_field(k0, force_amp);
		
		is_force_field_para_read = 1;  // To read only once.
	}	
	
}


//*********************************************************************************************
//
// Scalar
//

void IncFluid::Compute_force_Taylor_Green(IncSF& T)
{

	if (is_force_field_para_read == 0)
	{
		*T.Force = 0.0;
		
		Compute_force_Taylor_Green();
	
		is_force_field_para_read = 1;
	}	

}


//*********************************************************************************************
//
// force both V and W
//

void IncFluid::Compute_force_Taylor_Green(IncVF& W)
{
	
	if (is_force_field_para_read == 0)
	{
		int k0 = ((int) (*force_field_para)(1));
		DP force_amp = (*force_field_para)(2);
		DP forceW_amp = (*force_field_para)(3);
		
		Setup_Taylor_Green_force_field(k0, force_amp);
		W.Setup_Taylor_Green_force_field(k0, forceW_amp);
		
		is_force_field_para_read = 1;
	}
		
}


//*********************************************************************************************	
//
//  Force both V and W
//	
	
void IncFluid::Compute_force_Taylor_Green(IncVF& W, IncSF& T)
{

	if (is_force_field_para_read == 0)
	{
		*T.Force = 0.0;
		Compute_force_Taylor_Green(W);
		
		is_force_field_para_read = 1;
	}	
	
}	


//****************************** End of compute_force_TG.cc ***********************************



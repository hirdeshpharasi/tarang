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

// omega = 2 omega(orig) L/nu
void IncFluid::Compute_force_Coriolis()
{
	int omega_components = (*force_int_para)(1);
	
	if (omega_components == 1) {
		DP omega =  (*force_double_para)(1);
		
		if ((*force_int_para)(2) == 1) {	// omega along x
			*Force1 = 0.0;
			*Force2 = -omega*(*V3);
			*Force3 = omega*(*V2);
		}
		else if ((*force_int_para)(2) == 2) {// omega along y
			*Force1 = omega*(*V3);
			*Force2 = 0.0;
			*Force3 = -omega*(*V1);
		}
		else if ((*force_int_para)(2) == 3) {// omega along z
			*Force1 = -omega*(*V2);
			*Force2 = omega*(*V1);
			*Force3 = 0.0;
		}
	}
	
	else if (omega_components == 3) {
		DP omega1 = (*force_double_para)(1);
		DP omega2 = (*force_double_para)(2);
		DP omega3 = (*force_double_para)(3);
		
		*Force1 = omega2*(*V3) - omega3*(*V2);
		*Force2 = omega3*(*V1) - omega1*(*V3);
		*Force3 = omega1*(*V2) - omega2*(*V1);
	}	
}


//*********************************************************************************************
//
// Scalar
//

void IncFluid::Compute_force_Coriolis(IncSF& T)
{
	Compute_force_Coriolis();

}


//*********************************************************************************************
//
// force both V and W
//

void IncFluid::Compute_force_Coriolis(IncVF& W)
{
	Compute_force_Coriolis();
		
}


//*********************************************************************************************	
//
//  Force both V and W
//	
	
void IncFluid::Compute_force_Coriolis(IncVF& W, IncSF& T)
{
	Compute_force_Coriolis();
	
}	


//****************************** End of compute_force_TG.cc ***********************************



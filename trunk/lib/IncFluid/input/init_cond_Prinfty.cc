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


/*! \file  init_cond_Prinfty.cc
 * 
 * @brief Initial conditions for Pr=infty convection.
 *
 * @note Given T(k), find u(k).
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "../IncFluid.h"



//*********************************************************************************************


void IncFluid::Init_cond_Prinfty(IncSF& T)
{
	
	Compute_force_RB(T);  // F = R\theta \hat{x}
	
	Add_force(T);
	
	Compute_pressure();   // *F contains pressure
	
	Compute_rhs(T);
		
	*V1 = *nlin1;
	*V2 = *nlin2;
	*V3 = *nlin3;
	
	Array_divide_ksqr_SCFT(N, *V1, kfactor);				// u(k) = nlin(k)/k^2	
	Array_divide_ksqr_SCFT(N, *V2, kfactor);
	Array_divide_ksqr_SCFT(N, *V3, kfactor); 
}

//******************************** End of Init_cond_Prinfty.cc *************************************

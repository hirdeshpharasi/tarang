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


/*! \file  init_cond_TG.cc
 * 
 * @brief Initial conditions as Taylor Green flow (TG).
 *
 * @note The parameters are read from parameter file. 
 *
 *		Vx = amp*sin(k0 x) cos(k0 y) cos(k0 z)
 *		Vy = -amp*cos(k0 x)sin(k0 y) cos(k0 z)
 *		Vz = 0
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "../IncFluid.h"


//*********************************************************************************************


void IncFluid::Init_cond_Taylor_Green()
{
	
	int k0 = (*init_cond_int_para)(1);
	DP amp = (*init_cond_double_para)(1);
	
	Setup_Taylor_Green_field(k0, amp);
}

//
// Scalar 
//


void IncFluid::Init_cond_Taylor_Green(IncSF& T)
{
	if ((globalvar_prog_kind == "INC_SCALAR") || (globalvar_prog_kind == "INC_SCALAR_DIAG"))
		Init_cond_Taylor_Green_scalar(T);
	
	else if ((globalvar_prog_kind == "RB_SLIP") || (globalvar_prog_kind == "RB_SLIP_DIAG"))
		Init_cond_Taylor_Green_RB(T);
}


void IncFluid::Init_cond_Taylor_Green_scalar(IncSF& T)
{	
	Init_cond_Taylor_Green();
	
	*T.F = 0.0;
}


void  IncFluid::Init_cond_Taylor_Green_RB(IncSF& T)
{
	
	if (globalvar_Pr_switch == "PRZERO") 
	{
		Init_cond_Taylor_Green();	
		
		*T.F = *V1; 
		Array_divide_ksqr(basis_type, N, *T.F, kfactor);		
	}
	
	else if (globalvar_Pr_switch == "PRINFTY") 
	{
		cout << "TG forcing not allowed for PRINFTY option" << endl;
		exit(0);
	}
	
	else
		Init_cond_Taylor_Green(T);
	
	Zero_modes_RB_slip(T);			
}

//
// Vector
//

void IncFluid::Init_cond_Taylor_Green(IncVF& W)
{

	int k0 = (*init_cond_int_para)(1);
	DP amp = (*init_cond_double_para)(1);
	DP ampW = (*init_cond_double_para)(3);
	
	Setup_Taylor_Green_field(k0, amp);
	W.Setup_Taylor_Green_field_mag_field(k0, ampW);		
}

//
//	Vector+scalar
//

void IncFluid::Init_cond_Taylor_Green(IncVF& W, IncSF& T)
{	
	Init_cond_Taylor_Green(W);
	
	*T.F = 0.0;

}



//******************************** End of Init_cond_TG.cc *************************************








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

/*! \file  init_cond_RB_misc.cc
 * 
 * @brief Initial conditions as RB convection: Lorenz equations.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FOUR basis; not for SCFT basis.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 */

#include "../IncFluid.h"


//**********************************************************************************************


/** @brief Set the modes for Rayleigh Taylor instability
 * 
 * @return  rho = 0 for x=0:0.5, and 1 for x=0.5:1.
 *		in T:  1 for x=0:0.5, and 0 for x=0.5:1.
 *		we model by T(x) = ( 1-tanh((x-0.5)*slope) )/2 with slope = large no.
 *		theta = T(x) -1+x.
 */
void  IncFluid::Init_cond_Rayleigh_Taylor(IncSF& T)
{

	if (basis_type == "SCFT")
	{
		DP		slope = 1000;
		DP		x;
		
		*V1 = 0.0;
		*V2 = 0.0;
		*V3 = 0.0;
		
		for (int lx=0; lx<local_N1; lx++)
		{
			x = lx*1.0/(local_N1_start + lx);
			(*T.Fr)(lx,Range::all(),Range::all()) = x - (1+tanh(slope*(x-0.5)))/2;
		}

#ifndef TRANSPOSE
		T.RS_Forward_transform(*VF_temp_r);

#else
		T.RS_Forward_transform_transpose_order(*T.F, *VF_temp_r);	
		// Forward_tr(*T.Fr) = *T.F
#endif	
		
		
	}
	
	else
	{
		if (my_id == master_id)
		{
			cout << "Init_cond_Rayleigh_Taylor allowed not allowed for the chosen basis type" 
					<< endl;
			
			exit(1);
		}
	}
}


//******************************** End of Init_cond_RB_misc.cc ********************************




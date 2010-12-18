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

/*! \file  compute_force_main.cc
 * 
 * @brief  Compute force and put it in F.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"

//*********************************************************************************************


void IncFluid::Compute_force()
{
	switch (force_field_proc) 
	{		
		case (0) : Compute_force_decay();	break;				
		case (1) : Compute_force_const_energy_helicity_supply(); break;	
		case (2) : Compute_force_const_energy_helicity(); break;			
		case (3) : Compute_force_Taylor_Green(); break;
		case (4) : Compute_force_ABC(); break;
		case (5) : Compute_force_given_modes_SIMPLE(); break;
		case (6) : Compute_force_given_modes_VORTICITY(); break;		
		case (11): Compute_force_Liquid_metal(); break;
	}	
}

//
//

void IncFluid::Compute_force(IncSF& T)
{
	switch (force_field_proc) 
	{		
		case (0) : Compute_force_decay(T);	break;					
		case (1) : Compute_force_const_energy_helicity_supply(T); break;	
		case (2) : Compute_force_const_energy_helicity(T); break;		
		case (3) : Compute_force_Taylor_Green(T); break;
		case (4) : Compute_force_ABC(T); break;
		case (5) : Compute_force_given_modes_SIMPLE(T); break;
		case (6) : Compute_force_given_modes_VORTICITY(T); break;
		case (51) : Compute_force_RB(T); break;	
		case (52): Compute_force_nonboussenesq(T); break;	
	}
}

//
//

void IncFluid::Compute_force(IncVF& W)
{
	switch (force_field_proc) 
	{		
		case (0) : Compute_force_decay(W);	break;				
		case (1) : Compute_force_const_energy_helicity_supply(W); break;	
		case (2) : Compute_force_const_energy_helicity(W); break;		
		case (3) : Compute_force_Taylor_Green(W); break;
		case (4) : Compute_force_ABC(W); break;
		case (5) : Compute_force_given_modes_SIMPLE(W); break;
		case (6) : Compute_force_given_modes_VORTICITY(W); break;
		case (11): Compute_force_DYNAMO_SIX_MODE(W); break;
	}	
}

//
//

void IncFluid::Compute_force(IncVF& W, IncSF& T)
{
	switch (force_field_proc) 
	{		
		case (0) : Compute_force_decay(W, T);	break;					
		case (1) : Compute_force_const_energy_helicity_supply(W, T); break;	
		case (2) : Compute_force_const_energy_helicity(W, T); break;			
		case (3) : Compute_force_Taylor_Green(W, T); break;
		case (4) : Compute_force_ABC(W, T); break;
		case (5) : Compute_force_given_modes_SIMPLE(W, T); break;
		case (6) : Compute_force_given_modes_VORTICITY(W, T); break;
		case (51) : Compute_force_RB(W,T); break;
	}
}


//*******************************  End of compute_force_main.cc *******************************




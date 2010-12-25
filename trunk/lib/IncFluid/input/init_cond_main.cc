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

/*! \file  init_cond_main.cc
 * 
 * @brief  Initial condition main file.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug No known bugs
 */

#include "../IncFluid.h"

extern Uniform<DP> SPECrand;

// const int IC_MAX_NO_WAVENO = 40;


/**********************************************************************************************

		Sets up initial conditions of the Vector and Scalar field

***********************************************************************************************/


void IncFluid::Read_init_cond()
{
	switch (field_input_proc) 
	{
		
		case (1) : Init_cond();	break;					
		// read from field_in_file
		
		case (2) : Init_cond_reduced();	break;	
		// read from field_in_file with Nreduced D
		
		case (3) : Init_cond_modes_SIMPLE();	break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
		case (4) : Init_cond_modes_VORTICITY(); break;		
		// Modes - ki, Vx,  vorticity, Theta
		
		case (5) : Init_cond_energy_spectrum(); break;			
		// given energy spectrum
		
		case (6) : Init_cond_energy_helicity_spectrum(); break;	
		// given energy and hel spectrum
		
		case (7) : Init_cond_Taylor_Green(); break;
		
		case (8) : Init_cond_ABC(); break;
	
	}
	
	if ((free_slip_verticalwall_switch == 1) && (basis_type == "SCFT"))
		free_slip_verticalwall_field();
}

//
//

void IncFluid::Read_init_cond(IncSF& T)
{
	
	switch (field_input_proc) 
	{
		
		case (1) : Init_cond(T);	break;			
		// read from field_in_file
		
		case (2) : Init_cond_reduced(T);	break;		
		// read from field_in_file with Nreduced D
		
		case (3) : Init_cond_modes_SIMPLE(T);	break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
		case (4) : Init_cond_modes_VORTICITY(T); break;		
		// Modes - ki, Vx,  vorticity, Theta
		
		case (5) : Init_cond_energy_spectrum(T); break;				
		// given energy spectrum
		
		case (6) : Init_cond_energy_helicity_spectrum(T); break;	
		// given energy and hel spectrum
		
		case (7) : Init_cond_Taylor_Green(T); break;	
		// initialize only V field
		
		case (8) : Init_cond_ABC(T); break;			
		// initialize only V field
			
		case (51) : Init_cond_RB_Lorenz(T); break;
		// Lorenz initcond for RB convection
			
		case (511) : Init_cond_Rayleigh_Taylor(T); break;			
		// initialize For Rayleigh Taylor instability.	
		
	}
	
	if ((free_slip_verticalwall_switch == 1) && (basis_type == "SCFT"))
		free_slip_verticalwall_field(T);
}

//
//

void IncFluid::Read_init_cond(IncVF& W)
{
	switch (field_input_proc)
	{
		
		case (1) : Init_cond(W);	break;			
		// read from field_in_file
		
		case (2) : Init_cond_reduced(W); break;	
		// read from field_in_file with Nreduced D
		
		case (3) : Init_cond_modes_SIMPLE(W);	break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
		case (4) : Init_cond_modes_VORTICITY(W); break;		
		// Modes - ki, Vx,  vorticity, Theta
		
		case (5) : Init_cond_energy_spectrum(W); break;				
		// given energy spectrum
		
		case (6) : Init_cond_energy_helicity_spectrum(W); break;	
		// given energy and hel spectrum
		
		case (7) : Init_cond_Taylor_Green(W); break;	
		// initialize V, W field
		
		case (8) : Init_cond_ABC(W); break;			
		// initialize V,W  field
		
		case (101) : Init_cond_DYNAMO_SIX_MODE(W); break;
	}
	
	if ((free_slip_verticalwall_switch == 1) && (basis_type == "SCFT"))
		free_slip_verticalwall_field(W);
}


//
//

void IncFluid::Read_init_cond(IncVF& W, IncSF& T)
{
	switch (field_input_proc) 
	{
		
		case (1) : Init_cond(W, T);	break;						
		// read from field_in_file
		
		case (2) : Init_cond_reduced(W, T);	break;			
		// read from field_in_file with Nreduced D
		
		case (3) : Init_cond_modes_SIMPLE(W, T);	break;		
		// Modes - ki, Vx, (Vy:3D),Theta
		
		case (4) : Init_cond_modes_VORTICITY(W, T); break;		
		// Modes - ki, Vx,  vorticity, Theta
		
		case (5) : Init_cond_energy_spectrum(W, T); break;			
		// given energy spectrum
		
		case (6) : Init_cond_energy_helicity_spectrum(W, T); break;	
		// given energy and hel spectrum
		
		case (7) : Init_cond_Taylor_Green(W, T); break;	
		// initialize V, W field
		
		case (8) : Init_cond_ABC(W, T); break;			
		// initialize V, W field
	
	}
	
	if ((free_slip_verticalwall_switch == 1) && (basis_type == "SCFT"))
		free_slip_verticalwall_field(W, T);
}



//******************************** End of Init_cond_main.cc ***********************************




  

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

/*! \file  IncSF.cc
 * 
 * @sa	IncSF.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug   No known bug
 */



#include "IncSF.h"


//*********************************************************************************************

IncSF::IncSF
(
	int NN[], 
	string string_switches[], 
	Array<int,1> switches, 
	DP *prog_kfactor, 
	DP diff_coefficient,  
	DP hyper_diff_coefficient, 
	Array<int,1> misc_output_para
): 
		CSF(NN, string_switches, switches, prog_kfactor, misc_output_para), 
		RSF(NN, string_switches) 
{

	ISF_basis_type				= string_switches[1];
	ISF_alias_switch			= string_switches[2];
		
	hyper_diffusion_switch		= switches(1);
	
	diffusion_coefficient		= diff_coefficient;
	hyper_diffusion_coefficient = hyper_diff_coefficient;
	
	for (int i=1; i<=3; i++) 
	{
		NIs[i]=NN[i]; 
		ISF_kfactor[i] = prog_kfactor[i];
	}  
	

	nlin = new Array<complx,3>(local_N1, NIs[2], NIs[3]/2+1);	// nlin contains FT[u.grad(T)]
	SF_temp = new Array<complx,3>(local_N1, NIs[2], NIs[3]/2+1); 
	Force = new Array<complx,3>(local_N1, NIs[2], NIs[3]/2+1); 
	
	(*nlin) = 0.0; 
	(*SF_temp) = 0.0;	
	(*Force) = 0.0; 
	
	shell_spectrum_force_SFk = new Array<double,1>(CS_shell_ek_size);  
	(*shell_spectrum_force_SFk)=0.0;
	
	// Anisotropic spectrum
	
	if (CS_anisotropic_ring_switch == 1)
	{
		ring_spectrum_force_SFk = new Array<double,2>(CS_ring_ek_size, 
														CS_no_sectors_spectrum+1);  
		(*ring_spectrum_force_SFk)=0.0;
	}
	
	if 	(CS_anisotropic_cylinder_switch == 1)	
	{		
		cylinder_ring_spectrum_force_SFk = new Array<double,2>(CS_cylinder_shell_ek_size, 
															CS_no_cylinder_slabs_spectrum+1);  
		(*cylinder_ring_spectrum_force_SFk)=0.0;
	}
}

//********************************  End of IncSF.cc *******************************************




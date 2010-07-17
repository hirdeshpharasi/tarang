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

/*! \file  IncVF.cc
 * 
 * @brief  Class declaration of IncVF, Incompressible Vector Field 
 *
 * @sa IncVF.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No know bugs
 */

#include "IncVF.h"


//*********************************************************************************************


		
IncVF::IncVF
(
	int NN[], 
	string string_switches[], 
	Array<int,1> switches, 
	DP *prog_kfactor, 
	DP diss_coefficient, 
	DP hyper_diss_coefficient, 
	Array<int,1> misc_output_para,
	Array<int,1> ET_parameters, 
	Array<DP,1>  ET_shell_radii_sector_array
): 
	
		CVF(NN, string_switches, switches, prog_kfactor, misc_output_para), 
		RVF(NN, string_switches),
		NLIN(NN, string_switches, switches, prog_kfactor, misc_output_para), 
		EnergyTr(NN, string_switches, switches, prog_kfactor, 
					ET_parameters, ET_shell_radii_sector_array)
{
	
	basis_type		= string_switches[1];
	alias_switch	= string_switches[2];
	
	hyper_dissipation_switch = switches(1);
	
	sector_spectrum_input_scheme	= misc_output_para(1);
	no_sectors_spectrum				= misc_output_para(2);
	no_cylinder_slabs_spectrum		= misc_output_para(3);
	
	for (int i=1; i<=3; i++)
	{
		N[i]=NN[i]; 
		kfactor[i] = prog_kfactor[i];
	}
	
	dissipation_coefficient = diss_coefficient;
	hyper_dissipation_coefficient = hyper_diss_coefficient;
	
	
	if (CV_anisotropic_ring_switch == 1)
	{
		// Angles of the sectors

		sector_angle_array_spectrum = new Array<DP,1>(no_sectors_spectrum+1);
		
		*sector_angle_array_spectrum = *CV_sector_angle_array_spectrum;
	}
	
	if 	(CV_anisotropic_cylinder_switch == 1)	
	{	
		// Cylinder kz's
		
		cylinder_kpll_array_spectrum = new Array<DP,1>(no_cylinder_slabs_spectrum+1);
		
		*cylinder_kpll_array_spectrum = *CV_cylinder_kpll_array_spectrum;
	}
	
	
	// Memory allocation..
	  
	Force1 = new Array<complx,3>(local_N1, N[2],N[3]/2+1);
	Force2 = new Array<complx,3>(local_N1, N[2],N[3]/2+1);
	Force3 = new Array<complx,3>(local_N1, N[2],N[3]/2+1);
	
	(*Force1) = 0.0; 
	(*Force2) = 0.0; 
	(*Force3) = 0.0;
	
	VF_temp = new Array<complx,3>(local_N1, N[2],N[3]/2+1); 
	VF_temp2 = new Array<complx,3>(local_N1, N[2],N[3]/2+1);
	
	if (N[2] > 1)
		VF_temp_r = new Array<complx,3>(local_N2, N[1],N[3]/2+1);
	
	else if (N[2] == 1)
		VF_temp_r = new Array<complx,3>(local_N2, N[2],N[1]/2);
	

	(*VF_temp) = 0.0; 
	(*VF_temp2) = 0.0;
	(*VF_temp_r) = 0.0;
	
	
	shell_ek_cross_V1T = new Array<double,1>(CV_shell_ek_size);  
	shell_ek_cross_V2T = new Array<double,1>(CV_shell_ek_size);
	shell_ek_cross_V3T = new Array<double,1>(CV_shell_ek_size);
	(*shell_ek_cross_V1T) = 0.0;
	(*shell_ek_cross_V2T) = 0.0;
	(*shell_ek_cross_V3T) = 0.0;

	if (CV_anisotropic_ring_switch == 1)
	{
		ring_ek_cross_V1T = new Array<double,2>(CV_ring_ek_size, no_sectors_spectrum+1);  
		ring_ek_cross_V2T = new Array<double,2>(CV_ring_ek_size, no_sectors_spectrum+1);
		ring_ek_cross_V3T = new Array<double,2>(CV_ring_ek_size, no_sectors_spectrum+1);
		
		(*ring_ek_cross_V1T)=0.0;
		(*ring_ek_cross_V2T)=0.0;
		(*ring_ek_cross_V3T)=0.0;
	}
	
	if 	(CV_anisotropic_cylinder_switch == 1)
	{		
		cylinder_ring_ek_cross_V1T = new Array<double,2>(CV_cylinder_shell_ek_size, 
													no_cylinder_slabs_spectrum+1); 
													
		cylinder_ring_ek_cross_V2T = new Array<double,2>(CV_cylinder_shell_ek_size, 
													no_cylinder_slabs_spectrum+1);
													
		cylinder_ring_ek_cross_V3T = new Array<double,2>(CV_cylinder_shell_ek_size, 
													no_cylinder_slabs_spectrum+1);
		(*cylinder_ring_ek_cross_V1T)=0.0;
		(*cylinder_ring_ek_cross_V2T)=0.0;
		(*cylinder_ring_ek_cross_V3T)=0.0;
	}
	
	shell_spectrum_force_Vk = new Array<double,1>(CV_shell_ek_size);  
	(*shell_spectrum_force_Vk)=0.0;
	
	shell_ek_cross = new Array<double,1>(CV_shell_ek_size);  
	(*shell_ek_cross)=0.0;
	
	
	if (CV_anisotropic_ring_switch == 1)
	{
		ring_spectrum_force_Vk = new Array<double,2>(CV_ring_ek_size, no_sectors_spectrum+1);  
		(*ring_spectrum_force_Vk)=0.0;
		
		ring_ek_cross = new Array<double,2>(CV_ring_ek_size, no_sectors_spectrum+1);
		*ring_ek_cross = 0.0;
	}
	
	//
	
	if 	(CV_anisotropic_cylinder_switch == 1)
	{				
		cylinder_ring_spectrum_force_Vk = new Array<double,2>(CV_cylinder_shell_ek_size, 
													no_cylinder_slabs_spectrum+1); 
		(*cylinder_ring_spectrum_force_Vk)=0.0;

		
		cylinder_ring_ek_cross = new Array<double,2>(CV_cylinder_shell_ek_size, 
													no_cylinder_slabs_spectrum+1);
		*cylinder_ring_ek_cross = 0.0;
	}
		
}

//*************************** End Class defs of  IncVF  ***************************************	





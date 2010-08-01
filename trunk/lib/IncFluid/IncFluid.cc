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

/*! \file  IncFluid.cc
 * 
 * @brief  Class constructor of IncFluid, Incompressible fluid field. 
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug No known bugs
 */
 
 
#include "IncFluid.h" 

 
//*********************************************************************************************

IncFluid:: IncFluid
(
	int N[], 
	string string_switches[], 
	Array<int,1> switches, 
	DP prog_kfactor[],
	DP diss_coefficient, 
	DP hyper_diss_coefficient,	
	Array<DP,1> time_para, 
	Array<DP,1> time_save_interval, 
	Array<int,1> misc_output_para,  
	Array<int,1> no_output_k_r, 
	Array<int,2> out_k_r_array, 
	Array<int,1> ET_parameters, 
	Array<DP,1> ET_shell_radii_sector_array, 
	Array<int,1> field_input_meta_para, 
	Array<DP,1> init_cond_parameters,
	Array<int,1> force_field_meta_para, 
	Array<DP,1> force_field_parameters
): 
		IncVF(N, string_switches, switches, prog_kfactor, 
				diss_coefficient, hyper_diss_coefficient,  misc_output_para,
				ET_parameters, ET_shell_radii_sector_array),
		Time(time_para, time_save_interval) 				
{


	data_dir_name				= string_switches[0];

	integ_scheme				= string_switches[3];
	nos_input_field_mode		= string_switches[4];
	nos_output_field_mode		= string_switches[5];

	low_dim_switch				= switches(2);
	les_switch					= switches(3);
	skpq_switch					= switches(8);
	output_field_r_switch		= switches(9);
	sincos_horizontal_2D_switch  = switches(10);
	
	// string_switches_all and *switches_all are useful in time_advance functions only.	
	for (int i=0; i<MAXSIZE_STRING_SWITCHES_ARRAY; i++)
		string_switches_all[i] = string_switches[i];
	
	switches_all = 	new Array<int,1>(MAXSIZE_SWITCHES_ARRAY);	
	for (int i=0; i<MAXSIZE_SWITCHES_ARRAY; i++)
		(*switches_all)(i) = switches(i);
	
	misc_output_para_all = 	new Array<int,1>(MAXSIZE_MISC_OUTPUT_PARA);		
	for (int i=0; i<MAXSIZE_MISC_OUTPUT_PARA; i++)
		(*misc_output_para_all)(i) = misc_output_para(i);
	// Contains details on energy transfer shell, rings etc.
	
///	

	for (int i=1; i<=3; i++)
		N_out_reduced[i] = misc_output_para(5+i);
	
	nos_output_waveno	= no_output_k_r(1);	
	output_k_array		= new Array<int,2>(nos_output_waveno+1, 4);
	
	nos_output_position	= no_output_k_r(2);
	output_position_array= new Array<int,2>(nos_output_position+1, 4);
	
	*output_k_array = 0;
	*output_position_array = 0;
	
	(*output_k_array)(Range(1, nos_output_waveno), Range::all()) 
		= out_k_r_array(Range(1, nos_output_waveno), Range::all());
			
	(*output_position_array)(Range(1, nos_output_position), Range::all()) 
		= out_k_r_array(Range(nos_output_waveno+1, nos_output_waveno+nos_output_position), 
								Range::all());
	
	if (my_id == master_id)		
	{
		cout << "Output_k_array(-, kx, ky, kz) " 
			<< (*output_k_array)(Range(1, nos_output_waveno), Range::all()) 
			<< endl << endl;
			
		cout << "Output_position_array(-, ix, iy, iz) " 
			<< (*output_position_array)(Range(1, nos_output_position), Range::all()) 
			<< endl << endl;	
	}

	// Initial condition parameters
	
	field_input_proc			= field_input_meta_para(1);
	number_of_init_cond_para	= field_input_meta_para(2);
	
	
	for (int i=1; i<=3; i++)
		N_in_reduced[i] = field_input_meta_para(i+2);
	
	init_cond_para  = new Array<DP,1>(MAXSIZE_INIT_COND_FIELD_PARA);
	*init_cond_para = 0.0;							// initialize to 0.0
	
	*init_cond_para = init_cond_parameters;			// copy 
		
	
	// field_modes for initial condition
	field_modes_k_array = new Array<int,2>(INIT_COND_MAX_NO_FIELD_MODES,4);
	field_modes_V_array = new Array<complx,2>(INIT_COND_MAX_NO_FIELD_MODES, 
													MAX_NO_VARIABLES+1);
	(*field_modes_k_array) = 0; 
	(*field_modes_V_array) = 0.0;
	
	num_field_waveno = 0;							
	// gets nonzero value from the field_in file
	
	
	// Force field parameters.
	force_field_proc		= force_field_meta_para(1);
	number_of_force_para	=  force_field_meta_para(2);
	
	force_field_para = new Array<DP, 1>(MAXSIZE_FORCE_PARA);
	*force_field_para = 0.0;						// initialize to 0.0
	
	*force_field_para = force_field_parameters;		// Copy
	
	// force_field_modes to be read
	
	force_field_modes_k_array = new Array<int,2>(MAX_NO_FORCE_FIELD_MODES,4);
	force_field_modes_F_array = new Array<complx,2>(MAX_NO_FORCE_FIELD_MODES, 
													MAX_NO_VARIABLES+1);	
	(*force_field_modes_k_array) = 0;
	(*force_field_modes_F_array) = 0.0;
	
	is_force_field_para_read = 0;
	
	num_force_waveno = 0;
	
	
	if (low_dim_switch == 1)
	{
		cout << " low-dimensional simulations not allowed for Parallel program " << endl;
		cout << " PROGRAM EXITING " << endl;
		exit(1);
	}

}


//******************************   End of IncFluid.cc  ****************************************


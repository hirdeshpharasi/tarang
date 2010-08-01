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


/*! \file  main.h
 * 
 *	@brief  Main file for Tarang (MPI). Calls various modules.
 *
 *	@author  M. K. Verma
 *	@version 4.0  MPI
 *	@date Sept 2008
 */

//********************************************************************************************* 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include "IncFluid.h"   
// Contains IncFlow declarations (which contains fields +fourier+sincosfour defs)

// declarations of global const variables

void Read_prog_para(ifstream& prog_para_file, string& prog_kind, 
						string& data_dir_name);

void Read_para
(
	ifstream& para_file, 
	int dim, 
	int number_of_fields,  
	int N[], 
	string string_switches[], 
	Array<int,1> switches,
	double diss_coefficient[], 
	double hyper_diss_coefficient[],
	Array<int,1> solver_meta_para,
	Array<int,1> solver_int_para,
	Array<DP,1> solver_double_para,
	string solver_string_para[],
	Array<DP,1> time_para, 
	Array<DP,1> time_save, 
	Array<int,1> misc_output_para,
	Array<int,1> ET_parameters, 
	Array<DP,1> ET_shell_radii_sector_array,
	Array<int,1> no_output_k_r, 
	Array<int,2> out_k_r_array,
	Array<int,1> field_input_meta_para, 
	Array<DP,1> init_cond_para,
	Array<int,1> force_field_meta_para, 
	Array<DP,1> force_field_para
);


void Read_diag_para
(
 ifstream& para_file, 
 int dim, 
 int number_of_fields,  
 int N[], 
 string string_switches[], 
 Array<int,1> switches,
 double diss_coefficient[], 
 double hyper_diss_coefficient[],
 Array<int,1> solver_meta_para,
 Array<int,1> solver_int_para,
 Array<DP,1> solver_double_para,
 string solver_string_para[],
 Array<int,1>  diagnostic_procedure,
 Array<DP,1> time_para, 
 Array<int,1> misc_output_para,
 Array<int,1> ET_parameters, 
 Array<DP,1> ET_shell_radii_sector_array,
 Array<int,1> no_output_k_r, 
 Array<int,2> out_k_r_array,
 Array<int,1> field_input_meta_para, 
 Array<DP,1> init_cond_para,
 Array<int,1> force_field_meta_para, 
 Array<DP,1> force_field_para
);


														
int Ifluid_main(string data_dir_name);

int Ifluid_diag_main(string data_dir_name);

int Iscalar_main(string data_dir_name);

int Iscalar_diag_main(string data_dir_name);

int IMHD_main(string data_dir_name);

int IMHD_diag_main(string data_dir_name);

int RB_slip_main(string data_dir_name);

int RB_slip_diag_main(string data_dir_name);

int RB_slipMHD_main(string data_dir_name);

	
//********************************** End of main.h ********************************************



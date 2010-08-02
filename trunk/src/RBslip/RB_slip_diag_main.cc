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

/*! \file RB_slip_main.cc 
 * 
 * @brief Main program executing RB convection under free_slip bounadary condition
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */  


#include "RB_slip_main.h"

// Global variables 
extern int my_id;								// My process id
extern int numprocs;							// No of processors
extern const int master_id;						// Id of master proc
extern ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
extern ptrdiff_t local_N2, local_N2_start;
extern MPI_Status status;
					
extern Uniform<DP> SPECrand;					// for random vars--  declared in main.cc


//********************************************************************************************* 
 
int RB_slip_diag_main(string data_dir_name)
{
	
	cout << "ENTERING RB_slip_main :   my_id " << my_id << endl;
	
	//**************** Variable declarations *******************
	
	// Variables connected to the field variables and integration
	
	int			N[4];
	int			no_of_fields = 2;				// U, B fields
	DP			diss_coefficients[2];
	DP			hyper_diss_coefficients[2];
	
	
	ifstream	para_file;					// field and input-output parameters stored here
	
	string		string_switches[MAXSIZE_STRING_SWITCHES_ARRAY];
	
	Array<int,1>	switches( MAXSIZE_SWITCHES_ARRAY);
	
	// solver para
	Array<int,1>	solver_meta_para(MAXSIZE_SOLVER_META_PARA);
	Array<int,1>	solver_int_para(MAXSIZE_SOLVER_INT_PARA);
	Array<DP,1>		solver_double_para(MAXSIZE_SOLVER_DOUBLE_PARA);
	string			solver_string_para[MAXSIZE_SOLVER_STRING_PARA];

	
	Array<int,1>	diagnostic_procedure(MAXSIZE_DIAGNOSTIC_PARA);
	diagnostic_procedure = 0;
	
	// Variables connected to the output functions
	
	//! Tinit, Tfinal, Tdt, Tdiagnostic_int in index 1..4
	Array<DP, 1>time_para(MAXSIZE_TIME_PARA_ARRAY);					
	time_para = 0.0;

	Array<DP,1>		time_save(MAXSIZE_TIME_SAVE_ARRAY);
	time_save		= 0.0;
	
	Array<int,1>	misc_output_para(MAXSIZE_MISC_OUTPUT_PARA);
	
	//!  contains no_output_k, no_output_r
	Array<int,1>	no_output_k_r(3);							
	Array<int, 2>	out_k_r_array(MAXSIZE_OUT_K_R_ARRAY,4);
	out_k_r_array		= 0;			// intialize
	
	
	
	// Variables connected to the energy transfer
	
	Array<int,1>	ET_parameters(MAXSIZE_ET_PARAMETERS);
	
	Array<DP,1>		ET_shell_radii_sector_array(MAXSIZE_ET_RADII_SECTOR_ARRAY);
	
	ET_shell_radii_sector_array = 0.0;



	// Variables connected to the initial conditions
	
	Array<int,1>	field_input_meta_para(MAXSIZE_INIT_COND_FIELD_META_PARA);			
	
	Array<DP,1>		init_cond_para(MAXSIZE_INIT_COND_FIELD_PARA);	
	
	
	// Variables connected to the force field
		
	Array<int,1>	force_field_meta_para(MAXSIZE_FORCE_META_PARA);
	
	Array<DP, 1>	force_field_para(MAXSIZE_FORCE_PARA);
	
	
	// Variables connected to RB convection parameters

	DP			Pr;							// Prandtl number
	DP			Ra;							// Rayleigh number
	DP			r;							// r=Ra/Rac
	DP			k0;							// kfactor[2]=pi/sqrt(2)
	DP			qbyk0;						// q/k0
	
	string		Pr_switch;					// PRLARGE, PRSMALL, PRZERO
	string		RB_Uscaling;				// ULARGE, USMALL
		
	ifstream	RB_slip_para_file;			// RB para read from this file: data_dir_name/in/RB_para.d
	
	
	
	
	//**************** Program starts here ****************
	
	// Read field para and input-output parameters
	
	string		filename = "/in/para.d"; 
	filename = data_dir_name+ filename;  
	para_file.open(filename.c_str());		// filename = data_dir_name/in/para.d
	
	Read_diag_para(para_file, 3, no_of_fields, N, string_switches, switches,
				   diss_coefficients, hyper_diss_coefficients, 
				   solver_meta_para, solver_int_para, solver_double_para, solver_string_para,
				   diagnostic_procedure,
				   time_para,  
				   misc_output_para, ET_parameters, ET_shell_radii_sector_array,
				   no_output_k_r, out_k_r_array,
				   field_input_meta_para, init_cond_para,
				   force_field_meta_para, force_field_para);

	string_switches[0] = data_dir_name;
	
	globalvar_anisotropy_switch = switches(13);
	globalvar_waveno_switch = switches(14);
		
	
	// Read RB parameters like Pr, r, etc.
	
	filename = "/in/RB_para.d"; 
	filename = data_dir_name+ filename;					// filename = data_dir_name/in/RB_para.d	
	RB_slip_para_file.open(filename.c_str());
	
	Read_RB_slip_para(RB_slip_para_file, Pr, r, qbyk0,  Pr_switch, RB_Uscaling);

		
	// Construct output prefix for all the output files
	string prefix_str,  Pr_str, r_str;							
	ostringstream Pr_buffer, r_buffer;
	
	Pr_buffer << Pr;  
	Pr_str = Pr_buffer.str(); 
	
	r_buffer << r;	
	r_str = r_buffer.str();
	
	prefix_str = "%% Pr = " + Pr_str +  "  r = " + r_str;	// goes into all the output files		
	

	// kfactor computation
	Ra = (27.0*pow4(M_PI)/4.0)*r; 
	k0 = M_PI/sqrt(2.0);
	
	DP kfactor[4];
	
	/*	OLD SCHEME
	 kfactor[1] = M_PI; 
	kfactor[2] = k0; 
	q = kfactor[3] = k0 * qbyk0;
	 */
	
	
	if (solver_int_para(1) == 0)		// kfactor as given
	{
		kfactor[1] = solver_double_para(1);
		kfactor[2] = solver_double_para(2);
		kfactor[3] = solver_double_para(3);
			
	}
	
	else if (solver_int_para(1) == 1)	// kfactor multiplied by some predefined constant
	{
		kfactor[1] = M_PI * solver_double_para(1);
		kfactor[2] = k0 * solver_double_para(2);
		kfactor[3] = k0 * solver_double_para(3);
	}
	
	cout << "kfactor: " << kfactor[1] << " "  << kfactor[2] << " "  << kfactor[3] 
						<< endl << endl;
	
	
	// Set Dissipation coefficients
	
	if (Pr_switch == "PRZERO") 
	{
		diss_coefficients[0] = 1.0;
		diss_coefficients[1] = 0.0;
	}
	
	else if (Pr_switch == "PRLARGE") 
	{
		if (RB_Uscaling == "USMALL") 
		{
			diss_coefficients[0] = Pr;              //  Coeff of grad^2 u
			diss_coefficients[1]  = 1.0;			// Coeff of grad^2 T
		}
		else if (RB_Uscaling == "ULARGE") 
		{
			diss_coefficients[0] = sqrt(Pr/Ra);             
			diss_coefficients[1]  = 1/sqrt(Pr*Ra);			
		}
	}
	
	else if (Pr_switch == "PRSMALL") 
	{
		if (RB_Uscaling == "USMALL") 
		{
			diss_coefficients[0] = 1.0;             
			diss_coefficients[1]  = 1/Pr;			
		}
		else if (RB_Uscaling == "ULARGE") 
		{
			diss_coefficients[0] = sqrt(Pr/Ra);             
			diss_coefficients[1]  = 1/sqrt(Pr*Ra);			
		}
	}
	
	
	string basis_type  = string_switches[1];
	
	// Local_N1, local_N2 assignments 
	if (N[2] > 1)
	{
		if (basis_type == "FOUR")
		{
			
			int alloc_local;											  
			alloc_local = fftw_mpi_local_size_3d_transposed(N[1], N[2], N[3], MPI_COMM_WORLD,
										&local_N1, &local_N1_start, &local_N2, &local_N2_start);
		}
		
		else if (basis_type == "SCFT")
		{
			local_N1 = N[1]/numprocs;			// basis_type = SCFT
			local_N2 = N[2]/numprocs;
			local_N1_start = my_id * local_N1;
			local_N2_start = my_id * local_N2;		
		}
	}
	
	else if (N[2] == 1)
	{
		if (basis_type == "FOUR")
		{
			
			int alloc_local;											  
			alloc_local = fftw_mpi_local_size_2d_transposed(N[1], N[3], MPI_COMM_WORLD,
															&local_N1, &local_N1_start, &local_N2, &local_N2_start);
		}
		
		else if (basis_type == "SCFT")
		{
			local_N1 = N[1]/numprocs;			
			local_N2 = N[3]/numprocs;
			local_N1_start = my_id * local_N1;
			local_N2_start = my_id * local_N2;		
		}
	}	
	
	if (my_id == master_id) cout << "No or processors = " << numprocs << endl << endl;
	
	cout << "MY ID, local_N1, local_N1_start, local_N2, local_N2_start = " << my_id << "  "
		 << local_N1 << "  " << local_N1_start << "  " << local_N2 << "  " << local_N2_start
		 << endl << endl;


	   
	// Constructors of Vector and Scalar fields

	IncFluid  U(N, string_switches, switches, kfactor, 
					diss_coefficients[0], hyper_diss_coefficients[0], 
					time_para, time_save,
					misc_output_para, no_output_k_r, out_k_r_array,
					ET_parameters, ET_shell_radii_sector_array,
					field_input_meta_para, init_cond_para,
					force_field_meta_para, force_field_para
				);	
				
	IncSF     T(N, string_switches, switches, kfactor, 
					diss_coefficients[1], hyper_diss_coefficients[1],
					misc_output_para
				);
	
	
	// Create FFTW plans
	
	U.Init_fftw_plan();	
	
	// Initialize vector fields;
	
	if (my_id == master_id)		
		U.Open_field_input_files();
		
	U.Read_init_cond_RB(Pr_switch, T);	
	
	if (my_id == master_id)		
		U.Close_field_input_files();
			
	//
	//
	
	if (my_id == master_id)		
	{
		U.Open_output_files(prefix_str);
		cout << endl << "STARTING THE SIMULATION NOW" << endl;
	}


	/*
	 
	 Testting module..
	 ===============
	real(*U.V1) =  tensor::k+  (N[1]*(tensor::i));
	
	cout << "LOCAL N1, N3 " << local_N1 << " "  << local_N2 << endl << endl;
	
	cout << "stage 1 " << *U.V1 << " " << (*U.V1)(1,0,1) << " "  <<  (*U.V1)(0,0,1) << endl << endl;
	
	 Transpose_array(N,*U.V1, *U.VF_temp_r);
	
	cout << "stage 2 " << *U.VF_temp_r << " " << (*U.VF_temp_r)(1,0,1) << " "  <<  (*U.VF_temp_r)(0,0,1) << endl;
	
	Inverse_transpose_array(N, *U.VF_temp_r, *U.V1);

	cout << "stage 3 " << *U.VF_temp_r << " " << (*U.VF_temp_r)(1,0,1) << " "  <<  (*U.VF_temp_r)(0,0,1) << endl;
	

	cout << "stage 1 " << *U.V1 << endl << endl;  
	
	
		//	ArrayISFT_SCFT(isintr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, N, *U.V1, *U.VF_temp_r);
	
				 ArrayICFT_SCFT(icostr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, N, *U.V1, *U.VF_temp_r);
	
			cout << "After icft " << *U.V1 << " " << (*U.V1)(1,0,1) << " "  << (*U.V1)(2,0,1) << endl<< endl ;
	
		//	ArraySFT_SCFT(sintr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, N, *U.V1, *U.VF_temp_r);
	
		 ArrayCFT_SCFT(costr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, N, *U.V1, *U.VF_temp_r);
	
			cout << "After cft " << *U.V1 << " " << (*U.V1)(1,0,1) << " "  << (*U.V1)(2,0,1) << endl<< endl;
	
	
	exit(1);
	*/
	

	U.Tnow = U.Tinit;
	
	
	int i=1;
	while (diagnostic_procedure(i) < MY_MAX_INT)
	{
		switch (diagnostic_procedure(i)) 
		{		
			case (0) : U.Output_global(T, Ra, Pr, Pr_switch, RB_Uscaling);		break;					
			case (1) : U.Output_shell_spectrum(T, Pr_switch); break;	
			case (2) : U.Output_ring_spectrum(T, Pr_switch);  break;			
			case (3) : U.Output_cylinder_ring_spectrum(T, Pr_switch);  break;
			case (4) : U.Output_flux(T, Pr_switch);		break;
			case (5) : U.Output_shell_to_shell(T, Pr_switch); break;
			case (6) : U.Output_cylinder_ring_to_ring(T, Pr_switch);  break;		
			case (7) : U.Output_structure_fn(T, Pr_switch);  break;
			case (8) : U.Output_planar_structure_fn(T, Pr_switch);  break;
			case (9) : U.Output_moment(T, Pr_switch);  break;
			case (10) : U.Output_field_k(T, Pr_switch);  break;
			case (11) : U.Output_field_r(T, Pr_switch);  break;
			case (12) : U.Output_Skpq(T, Pr_switch);  break;
			case (13) : U.Output_realfield(T, Pr_switch);  break;
			case (14) : U.Output_field_reduced(T, Pr_switch); break;
		}	
		
		i++;
	}	
	

	if (my_id == master_id)		
		U.Close_output_files();
  
  
	return(1);
  
} 


//**************************** End of RB_slip_main.cc *****************************************




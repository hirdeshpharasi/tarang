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

/*! \file Ifluid_main.cc 
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author  M. K. Verma
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 


 #include "Ifluid_main.h"
 
 // Global variable 
extern int my_id;								// My process id
extern int numprocs;							// No of processors
extern const int master_id;						// Id of master proc
extern ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
extern ptrdiff_t local_N2, local_N2_start;
extern MPI_Status status;


extern Uniform<DP> SPECrand;					// for random vars--  declared in main.cc	

//********************************************************************************************* 					
 
int Ifluid_main(string data_dir_name)
{

	//**************** Variable declarations *******************
	
	// Variables connected to the field variables and integration
	
	cout << "ENTERING IFLUID MAIN my_id " << my_id << endl;
	
	int			N[4];		
	int			no_of_fields = 1;				// U field
	DP			diss_coefficients[1];
	DP			hyper_diss_coefficients[1];
	
	ifstream	para_file;					// field and input-output parameters stored here
	
	string		string_switches[MAXSIZE_STRING_SWITCHES_ARRAY];
	
	Array<int,1>	switches( MAXSIZE_SWITCHES_ARRAY);
	
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

	
	//**************** Program starts here ****************
		
	// Read field para, input-output, force field parameters
	
	string		filename = "/in/para.d"; 
	filename = data_dir_name+ filename;  
	cout << "file name " << filename << endl;
	
	para_file.open(filename.c_str());		// filename = data_dir_name/in/para.d
	
	Read_para(para_file, 3, no_of_fields, N, string_switches, switches,
					diss_coefficients, hyper_diss_coefficients, time_para, time_save, 
					misc_output_para, ET_parameters, ET_shell_radii_sector_array,
					no_output_k_r, out_k_r_array,
					field_input_meta_para, init_cond_para,
					force_field_meta_para, force_field_para);
					
	string_switches[0] = data_dir_name;
																					
											
	// Construct output prefix for all the output files
	string prefix_str,  nu_str;							
	ostringstream nu_buffer;
	
	nu_buffer << diss_coefficients[0];  
	nu_str = nu_buffer.str(); 	
	prefix_str = "%% nu = " + nu_str;					// goes into all the output files
	
	// kfactor computations	
	DP	kfactor[4];
	kfactor [1] = kfactor[2] = kfactor[3] = 1.0; 
	
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
			local_N1 = N[1]/numprocs;			
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
	
	if (my_id == master_id) 
		cout << "No or processors = " << numprocs << endl << endl;
	
	cout << "MY ID, local_N1, local_N1_start, local_N2, local_N2_start = " << my_id << "  "
		<< local_N1 << "  " << local_N1_start << "  " << local_N2 << "  " << local_N2_start
		<< endl << endl;
	
	// Constructors of Vector field

	IncFluid  U(N, string_switches, switches, kfactor, 
					diss_coefficients[0], hyper_diss_coefficients[0], 
					time_para, time_save,
					misc_output_para, no_output_k_r, out_k_r_array,
					ET_parameters, ET_shell_radii_sector_array,
					field_input_meta_para, init_cond_para,
					force_field_meta_para, force_field_para
				);
	
	// Create FFTW plans
																			
	 U.Init_fftw_plan();	
		
	// Initialize vector fields;
	
	if (my_id == master_id)		
		U.Open_field_input_files();	
		
	U.Read_init_cond();	
		
	if (my_id == master_id)	 
		U.Close_field_input_files();
	

	//
	//
	
	if (my_id == master_id)  
	{
		U.Open_output_files(prefix_str);
		cout << endl << "STARTING THE SIMULATION NOW" << endl;
	}


	int  iter=0;  // iterations 

	
	U.Tnow = U.Tinit;
	U.Output_all_inloop();  // for initial cond
	
		//	U.Open_output_files_hdf5(prefix_str, U.Tnow, N, 3);
		// U.Output_all_inloop_hdf5();

	do 
	{
		
		U.Compute_force();	

		U.Compute_nlin();
		U.Output_field_k_inloop();				// T(k) in the output computation needs nlin
												// Output at the present time.
													
		U.Satisfy_reality_condition();
		
		U.Add_force();		
		
		U.Compute_pressure();  
		U.Output_pressure_spectrum_inloop();		
		// Output pressure spectrum at the present time
		
		U.Tdt = U.Get_dt();
		U.Tnow = U.Tnow + U.Tdt;
		iter++;

		U.Time_advance();
		// FIELD AT new time.

		// Output at new time..
		U.CV_Compute_totalenergy_diss();
		if ( (my_id == 0) && isnan(U.CV_total_energy) ) 
		{ 
			cout << "ERROR: Numerical Overflow " << endl;  break; 
		}
		
		U.Output_all_inloop();

		//U.Open_output_files_hdf5(prefix_str, U.Tnow, N, 3);
		//U.Output_all_inloop_hdf5();
	}
	while (U.Tnow < U.Tfinal);
  
	if (my_id == master_id)	
		U.Close_output_files();
  
	return(1);
  
} 


//********************************** End of Ifluid_main.cc ************************************	



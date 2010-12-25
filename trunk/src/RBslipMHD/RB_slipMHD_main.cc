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


 #include "RB_slipMHD_main.h"
 
// Global variables 
extern int my_id;								// My process id
extern int numprocs;							// No of processors
extern const int master_id;						// Id of master proc
extern ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
extern ptrdiff_t local_N2, local_N2_start;
extern MPI_Status status;
					
extern Uniform<DP> SPECrand;					// for random vars--  declared in main.cc


//********************************************************************************************* 

 
int RB_slipMHD_main(string data_dir_name)
{
	cout << "ENTERING RB_slipMHD_main :   my_id " << my_id << endl;
	
	//**************** Variable declarations *******************
	
	// Variables connected to the field variables and integration
	
	int			N[4]; 
	int			no_of_fields = 3;			// U, B, T field
	DP			diss_coefficients[3];
	DP			hyper_diss_coefficients[3];
	
	ifstream	para_file;					// field and input-output parameters stored here
	
	string		string_switches[MAXSIZE_STRING_SWITCHES_ARRAY];
	
	Array<int,1>	switches( MAXSIZE_SWITCHES_ARRAY);
	
	
	// solver para
	Array<int,1>	solver_meta_para(MAXSIZE_SOLVER_META_PARA);
	Array<int,1>	solver_int_para(MAXSIZE_SOLVER_INT_PARA);
	Array<DP,1>		solver_double_para(MAXSIZE_SOLVER_DOUBLE_PARA);
	string			solver_string_para[MAXSIZE_SOLVER_STRING_PARA];
	
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
	DP			eta;						// magnetic diffusivity
	
	string		Pr_switch;					// PRLARGE, PRSMALL, PRZERO
	string		RB_Uscaling;				// ULARGE, USMALL
		
	ifstream RB_slip_para_file;				
	// RB para read from this file: data_dir_name/in/RB_para.d
	
	
	
	//**************** Program starts here ****************
	
	// Read field para and input-output parameters
	
	string		filename = "/in/para.d"; 
	filename = data_dir_name+ filename;  
	para_file.open(filename.c_str());		// filename = data_dir_name/in/para.d
	
	Read_para(para_file, 3, no_of_fields, N, string_switches, switches,
					diss_coefficients, hyper_diss_coefficients, 
					solver_meta_para, solver_int_para, solver_double_para, solver_string_para,
					time_para, time_save, 
					misc_output_para, ET_parameters, ET_shell_radii_sector_array,
					no_output_k_r, out_k_r_array,
					field_input_meta_para, init_cond_para,
					force_field_meta_para, force_field_para);

	string_switches[0] = data_dir_name;
	
	// Read RB parameters like Pr, r, etc.
	
	filename = "/in/RB_para.d"; 
	filename = data_dir_name+ filename;					// filename = data_dir_name/in/RB_para.d	
	RB_slip_para_file.open(filename.c_str());
	
	Read_RB_slipMHD_para(RB_slip_para_file, Pr, r, eta, qbyk0, Pr_switch, RB_Uscaling);

		
	// Construct output prefix for all the output files
	string prefix_str,  Pr_str, r_str, eta_str;							
	ostringstream Pr_buffer, r_buffer, eta_buffer;
	
	Pr_buffer << Pr;  
	Pr_str = Pr_buffer.str(); 
	
	r_buffer << r;	
	r_str = r_buffer.str();
	
	eta_buffer << eta;	
	eta_str = eta_buffer.str();
	
	prefix_str = "%% Pr = " + Pr_str +  "  r == " + r_str + "eta = " + eta_str;	
	// goes into all the output files
				
	Ra = (27.0*pow4(M_PI)/4.0)*r; 
	k0 = M_PI/sqrt(2.0);
	
	DP kfactor[4];  
	
	kfactor [1] = M_PI; 
	kfactor[2] = k0; 
	kfactor[3] = k0 * qbyk0;
	
	
	// Set Dissipation coefficients
	
	if (Pr_switch == "PRZERO") 
	{
		diss_coefficients[0] = 1.0;
		diss_coefficients[2] = 0.0;
	}
	
	else if (Pr_switch == "PRLARGE") 
	{
		if (RB_Uscaling == "USMALL") 
		{
			diss_coefficients[0] = Pr;              //  Coeff of grad^2 u
			diss_coefficients[2]  = 1.0;			// Coeff of grad^2 T
		}
		else if (RB_Uscaling == "ULARGE") 
		{
			diss_coefficients[0] = sqrt(Pr/Ra);             
			diss_coefficients[2]  = 1/sqrt(Pr*Ra);			
		}
	}
	
	else if (Pr_switch == "PRSMALL") 
	{
		if (RB_Uscaling == "USMALL") 
		{
			diss_coefficients[0] = 1.0;             
			diss_coefficients[2]  = 1/Pr;			
		}
		else if (RB_Uscaling == "ULARGE") 
		{
			diss_coefficients[0] = sqrt(Pr/Ra);             
			diss_coefficients[2]  = 1/sqrt(Pr*Ra);			
		}
	}

	diss_coefficients[1]  = eta;				// magnetic diffusivity.
	hyper_diss_coefficients[1] = 0.0;
	   
	// Constructors of Vector and Scalar fields

	IncFluid  U(N, string_switches, switches, kfactor, 
					diss_coefficients[0], hyper_diss_coefficients[0], 
					time_para, time_save,
					misc_output_para, no_output_k_r, out_k_r_array,
					ET_parameters, ET_shell_radii_sector_array,
					field_input_meta_para, init_cond_para,
					force_field_meta_para, force_field_para
				);	
	
	
	IncVF     W(N, string_switches, switches, kfactor,
					diss_coefficients[1], hyper_diss_coefficients[1],
					misc_output_para, 
					ET_parameters, ET_shell_radii_sector_array
				);
				
																
	IncSF     T(N, string_switches, switches, kfactor, 
					diss_coefficients[2], hyper_diss_coefficients[2],
					misc_output_para
				);
					
	
	// Create FFTW plans
	
	U.Init_fftw_plan();	
	
	// Initialize vector fields;
	
	if (my_id == master_id)		
		U.Open_field_input_files();
		
	U.Read_init_cond_RB(Pr_switch, W, T);	
	
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
	U.Output_all_inloop(W, T, Ra, Pr, Pr_switch, RB_Uscaling);		// for init cond
	
	do
	{
	
		U.Compute_force_RB(W, T, Ra, Pr, Pr_switch, RB_Uscaling);				
		// Both for V and Temperature field
		
		U.Compute_nlin(W, T, Pr_switch);		
		U.Output_field_k_inloop(W, T, Pr_switch);						
		// T(k) in the output computation needs nlin at the present time
		
		U.Add_force(W, T, Pr_switch);															
		
		U.Compute_pressure();  
		U.Output_pressure_spectrum_inloop();						
		// Output pressure spectrum at the present time
											
		U.Tdt = U.Get_dt(W, T);
		U.Tnow = U.Tnow + U.Tdt;
		iter++;
	
		U.Time_advance(W, T, Ra, Pr, Pr_switch, RB_Uscaling);
		// FIELD AT new time.
		
		// Output at new time..
		U.CV_Compute_totalenergy_diss(); 
		T.CS_Compute_totalenergy_diss();
		if ( isnan(U.CV_total_energy) || isnan(T.CS_total_energy) ) 
		{ 
			cout << "ERROR: Numerical Overflow " << endl;  
			break; 
		}
		
		if ((U.free_slip_verticalwall_switch == 1) && (U.basis_type == "SCFT"))
			U.free_slip_verticalwall_field(W, T);
		
		if (U.apply_realitycond_alltime_switch == 1)
			U.Satisfy_reality_condition_field(W. T);
		
		U.Output_all_inloop(W, T, Ra, Pr, Pr_switch, RB_Uscaling);
		
	}
	while (U.Tnow < U.Tfinal);
  
	if (my_id == master_id)
		U.Close_output_files();
  
} 




//**************************** End of RB_slipMHD_main.cc **************************************




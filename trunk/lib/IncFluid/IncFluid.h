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


/*! \file  IncFluid.h
 * 
 *	@brief  Class declaration of IncFluid, Incompressible fluid (V field).
 * 
 *	For fluid, scalar, MHD, this class is for velocity. <BR>
 *	For scalar simulation, we need IncFluid V & IncSF T. <BR>
 *	For MHD simulation, we need IncFluid V & IncVF W.
 *
 *	@author  M. K. Verma
 *	@version 4.0 MPI
 *	@date Sept 2008
 *
 * @bug  No known bugs
 */
 
//*********************************************************************************************

#include "../IncFlow/IncVF.h"   // Contains Incompressible scalar field IncVF defs
#include "IncFluid_Time.h"
#include "IncFluid_constants.h"

#ifndef _H_IncFluid
#define _H_IncFluid

/*
#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

*/
//! @brief Incompressible fluid
/*!
 *  Inherits IncVF, which itself inherits CVF, RVF, NLIN, and EnergyTr. <BR>
 *  Inherits Time that contains time vars. <BR>
 * 
 *	This is fluid 
 *
 *	@sa IncSF.h
 *	@sa Nlin.h
 *  @sa EnergyTr.h
 */
 
 
//*********************************************************************************************	

class IncFluid: public IncVF, public Time  // Vectorfield
{
public:

	
	string			data_dir_name;
	
	string			integ_scheme;
	string			nos_output_field_mode;			// ASCII or BINARY
	string			nos_input_field_mode;

/*
	string 			filename_hdf5;
	
	H5File* 		file;


	DataSet* 		realfield_dataset1;
	DataSet* 		realfield_dataset2;
	DataSet* 		realfield_dataset3;
	DataSet* 		field_dataset1R;
	DataSet* 		field_dataset1C;
	DataSet* 		field_dataset2R;
	DataSet* 		field_dataset2C;
	DataSet* 		field_dataset3R;
	DataSet* 		field_dataset3C;
*/
	int				low_dim_switch;					// 0: Not implemented for PLL
	int				les_switch;						// 1 if les on
	int				skpq_switch;
	int				output_field_r_switch;
	int				sincos_horizontal_2D_switch;	
	// 1 if horizontal basis function is cos(k0 z) or sin(k0 z).
	
	string			string_switches_all[MAXSIZE_STRING_SWITCHES_ARRAY];
	Array<int,1>	*switches_all;
	Array<int,1>	*misc_output_para_all;
	
	ifstream		field_in_file;
	ifstream		force_field_in_file;
	
	ofstream		global_file;
	ofstream		field_out_file;
	ofstream		field_frequent_out_file;
	ofstream		field_out_reduced_file;
	ofstream		realfield_out_file;
	ofstream		field_k_out_file;
	ofstream		field_r_out_file;
	ofstream		spectrum_file;
	ofstream		pressure_file;
	ofstream		flux_file;
	ofstream		shell_to_shell_file;		
	ofstream		ring_spectrum_file;
	ofstream		ring_to_ring_file;
	ofstream		cylinder_ring_spectrum_file;
	ofstream		cylinder_ring_to_ring_file;	
	ofstream		structure_fn_file;
	ofstream		planar_structure_fn_file;
	ofstream		moment_file;
	ofstream		Skpq_file;
	
	ofstream		readme_file;
	
	// Output
	int				N_out_reduced[4];
	
	// Output V(k) etc.; Max no of allowed wavenos is 99; 
	// contains 1-4:  kx, ky, kz
	int				nos_output_waveno;
	Array<int,2>	*output_k_array;					
	
														
	// Output V(x)
	int				nos_output_position;
	Array<int,2>	*output_position_array;					 
	
			

	// Initial condition parameters
	
	int				field_input_proc;
	int				number_of_init_cond_para;
	int				N_in_reduced[4];	
	Array<DP,1>		*init_cond_para;
		
	// field_modes for initial condition
	
	Array<int,2>	*field_modes_k_array;
	Array<complx,2> *field_modes_V_array;			// Contains the modes read
	
	int				num_field_waveno;				// No of vector field wavenumber read
													// as input of V.

	
	// Force parameters
	int				force_field_proc;				//!< forcing procedure, eg, 0 for decaying		
	int				number_of_force_para;
	Array<DP,1>		*force_field_para;
		
	// force_field_modes to be read
			
	Array<int,2>	*force_field_modes_k_array;	
	Array<complx,2> *force_field_modes_F_array;

	int				num_force_waveno;				// No of force_field wavenumbers
	int				is_force_field_para_read;		// 1 if read already; 0 otherwise.
											


//*********************************************************************************************		

public:
			
	IncFluid(int N[], string string_switches[], Array<int,1> switches, 
				DP prog_kfactor[],
				DP diss_coefficient, DP hyper_diss_coefficient,	
				Array<DP,1> time_para, Array<DP,1> time_save_interval, 
				Array<int,1> misc_output_para, 
				Array<int,1> no_output_k_r, Array<int,2> out_k_r_array, 
				Array<int,1> ET_parameters, Array<DP,1> ET_shell_radii_sector_array, 
				Array<int,1> field_input_meta_para, Array<DP,1> init_cond_parameters,
				Array<int,1> force_field_meta_para, Array<DP,1> force_field_parameters);
	
	void Compute_force();
	void Compute_force(IncSF& T);
	void Compute_force(IncVF& W);
	void Compute_force(IncVF& W, IncSF& T);		
					
	void Compute_force_RB(IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling);	
	void Compute_force_RB(IncVF& W, IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling);	
	
	void Compute_force_decay();
	void Compute_force_decay(IncSF& T);
	void Compute_force_decay(IncVF& W);
	void Compute_force_decay(IncVF& W, IncSF& T);		
	
	void Compute_force_const_energy_helicity_supply();
	void Compute_force_const_energy_helicity_supply(IncSF& T);
	void Compute_force_const_energy_helicity_supply(IncVF& W);
	void Compute_force_const_energy_helicity_supply(IncVF& W, IncSF& T);
	
	void Compute_force_const_energy_helicity();			
	void Compute_force_const_energy_helicity(IncSF& T);
	void Compute_force_const_energy_helicity(IncVF& W);						
	void Compute_force_const_energy_helicity(IncVF& W, IncSF& T);
	
	
	void Compute_force_Taylor_Green();
	void Compute_force_Taylor_Green(IncSF& T);
	void Compute_force_Taylor_Green(IncVF& W);
	void Compute_force_Taylor_Green(IncVF& W, IncSF& T);
	
	void Compute_force_ABC();
	void Compute_force_ABC(IncSF& T);
	void Compute_force_ABC(IncVF& W);
	void Compute_force_ABC(IncVF& W, IncSF& T);
	
	void Compute_force_DYNAMO_SIX_MODE(IncVF& W);
	void Compute_force_Liquid_metal();	
	
	void Read_waveno_force_field(int k_dim, int F_n, int& num_force_waveno);
	void Compute_force_given_modes_SIMPLE();
	void Compute_force_given_modes_SIMPLE(IncSF& T);
	void Compute_force_given_modes_SIMPLE(IncVF& W);
	void Compute_force_given_modes_SIMPLE(IncVF& W, IncSF& T);
	
	void Compute_force_given_modes_VORTICITY();
	void Compute_force_given_modes_VORTICITY(IncSF& T);
	void Compute_force_given_modes_VORTICITY(IncVF& W);
	void Compute_force_given_modes_VORTICITY(IncVF& W, IncSF& T);
	
			
	
	void Add_force();
	void Add_force(IncSF& T);
	void Add_force(IncVF& W);
	void Add_force(IncVF& W, IncSF& T);	
	void Add_force(IncSF& T, string Pr_switch);
	void Add_force(IncVF& W, IncSF& T, string Pr_switch);
	
	void Add_pressure_gradient();
			
	void Compute_rhs();       // nlin[i]=-ik[i]*pressure(k)+nlin[i] 
	void Compute_rhs(IncSF& T); 
	void Compute_rhs(IncVF& W);  // nlin[i]=-FT(grad(pressure)) +nlin[i]; W.nlin[i] = -W.nlin[i]
	void Compute_rhs(IncVF& W, IncSF& T); 
	void Compute_rhs(IncSF& T, string Pr_switch); 
	void Compute_rhs(IncVF& W, IncSF& T, string Pr_switch);
	
	DP Get_dt();
	DP Get_dt(IncSF& T);
	DP Get_dt(IncVF& W);
	DP Get_dt(IncVF& W, IncSF& T);
	
	
	void Single_time_step_EULER(DP dt);           
	void Single_time_step_EULER(IncSF& T, DP dt); 
	void Single_time_step_EULER(IncVF& W, DP dt); 
	void Single_time_step_EULER(IncVF& W, IncSF& T,  DP dt); 
	void Single_time_step_EULER(IncSF& T, string Pr_switch, DP dt);
	void Single_time_step_EULER(IncVF& W, IncSF& T, string Pr_switch, DP dt);
  
	void Single_time_step_RK2(DP dt);           // Time-advace V by a single-time-step
	void Single_time_step_RK2(IncSF& T, DP dt); // Time-advance V and Scalar T by a single-time-step
	void Single_time_step_RK2(IncVF& W, DP dt); // Time-advance V and Vector W by a single-time-step
	void Single_time_step_RK2(IncVF& W, IncSF& T,  DP dt); 
	void Single_time_step_RK2(IncSF& T, string Pr_switch, DP dt);
	void Single_time_step_RK2(IncVF& W, IncSF& T, string Pr_switch, DP dt);
	
	void Single_time_step_Semi_implicit(DP dt);           
	void Single_time_step_Semi_implicit(IncSF& T, DP dt); 
	void Single_time_step_Semi_implicit(IncVF& W, DP dt); 
	void Single_time_step_Semi_implicit(IncVF& W, IncSF& T,  DP dt); 
	void Single_time_step_Semi_implicit(IncSF& T, string Pr_switch, DP dt);
	void Single_time_step_Semi_implicit(IncVF& W, IncSF& T, string Pr_switch, DP dt);

	  
	void Time_advance();			// time advance 
	void Time_advance(IncSF& T);	//  time advance both V and scalar field T 
	void Time_advance(IncVF& W);	//  time advance U and B (MHD) 
  
  
	void Time_advance(IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling);
	void Time_advance(IncVF& W, IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling);
	
	/// File operations
	
	void Open_field_input_files();
	void Open_force_waveno_input_files();
	void Close_field_input_files();
	void Close_force_waveno_input_files();
	
	void Input_prefix(ifstream& field_in_file);
	void Input_force_prefix(ifstream& force_field_in_file);
	void Output_prefix(ofstream& X, string prefix_str);
	
	void Open_output_files(string prefix_str);
	void Open_output_files_hdf5(string prefix_str, float timeNow, int N[], int rank);
	void Close_output_files();
	void Close_output_files_hdf5(string key);
	
	// Initial condition
	
	void Read_init_cond();									
	void Read_init_cond(IncSF& T);	
	void Read_init_cond(IncVF& W);
	void Read_init_cond(IncVF& W, IncSF& T);
	void Read_init_cond_RB(string Pr_switch, IncSF& T);
	void Read_init_cond_RB(string Pr_switch, IncVF& W, IncSF& T);
	
		
	void Init_cond();
	void Init_cond(IncSF& T);
	void Init_cond(IncVF& W);
	void Init_cond(IncVF& W, IncSF& T);
	
	void Init_cond_reduced();
	void Init_cond_reduced(IncSF& T);
	void Init_cond_reduced(IncVF& W);
	void Init_cond_reduced(IncVF& W, IncSF& T);
	
	void Read_waveno_field(int k_dim, int V_n, int& num_field_waveno);
	void Init_cond_modes_SIMPLE();
	void Init_cond_modes_SIMPLE(IncSF& T);
	void Init_cond_modes_SIMPLE(IncVF& W);
	void Init_cond_modes_SIMPLE(IncVF& W, IncSF& T);
	
	void Init_cond_modes_VORTICITY();
	void Init_cond_modes_VORTICITY(IncSF& T);
	void Init_cond_modes_VORTICITY(IncVF& W);
	void Init_cond_modes_VORTICITY(IncVF& W, IncSF& T);
	
	void Init_cond_energy_spectrum();
	void Init_cond_energy_spectrum(IncSF& T);
	void Init_cond_energy_spectrum(IncVF& W);
	void Init_cond_energy_spectrum(IncVF& W, IncSF& T);
	
	
	void Init_cond_energy_helicity_spectrum();
	void Init_cond_energy_helicity_spectrum(IncSF& T);
	void Init_cond_energy_helicity_spectrum(IncVF& W);
	void Init_cond_energy_helicity_spectrum(IncVF& W, IncSF& T);
	
	
	void Init_cond_Taylor_Green();
	void Init_cond_Taylor_Green(IncSF& T);
	void Init_cond_Taylor_Green(IncVF& W);
	void Init_cond_Taylor_Green(IncVF& W, IncSF& T);
	
	void Init_cond_ABC();
	void Init_cond_ABC(IncSF& T);
	void Init_cond_ABC(IncVF& W);
	void Init_cond_ABC(IncVF& W, IncSF& T);

	void Init_cond_DYNAMO_SIX_MODE(IncVF& W);
	
	// Init condition RB Convection
	void Init_cond_RB_Lorenz(IncSF& T);  // Lorenz
	
	void Init_cond(string Pr_switch, IncSF& T);
	void Init_cond_reduced(string Pr_switch, IncSF& T);
	void Init_cond_modes_SIMPLE(string Pr_switch, IncSF& T);
	void Init_cond_modes_VORTICITY(string Pr_switch, IncSF& T);
	void Init_cond_energy_spectrum(string Pr_switch, IncSF& T);
	void Init_cond_energy_helicity_spectrum(string Pr_switch, IncSF& T);
	void Init_cond_Taylor_Green(string Pr_switch, IncSF& T);
	void Init_cond_ABC(string Pr_switch, IncSF& T);
	
	void Init_cond(string Pr_switch, IncVF& W, IncSF& T);
	void Init_cond_reduced(string Pr_switch, IncVF& W, IncSF& T);
	void Init_cond_modes_SIMPLE(string Pr_switch, IncVF& W, IncSF& T);
	void Init_cond_modes_VORTICITY(string Pr_switch, IncVF& W, IncSF& T);
	void Init_cond_energy_spectrum(string Pr_switch, IncVF& W, IncSF& T);
	void Init_cond_energy_helicity_spectrum(string Pr_switch, IncVF& W, IncSF& T);
	void Init_cond_Taylor_Green(string Pr_switch, IncVF& W, IncSF& T);
	void Init_cond_ABC(string Pr_switch, IncVF& W, IncSF& T);
	
	
  // Output functions
  
	void Output_all_inloop();
	void Output_all_inloop(IncSF& T);
	void Output_all_inloop(IncVF& W);
	void Output_all_inloop(IncVF& W, IncSF& T);
	
	void Output_all_inloop(IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling);
	void Output_all_inloop(IncVF& W, IncSF& T,  DP Ra, DP Pr, string Pr_switch, 
									string RB_Uscaling);
	void Output_all_inloop_hdf5();

	void Output_field_k_inloop();
	void Output_field_k_inloop(IncSF& T);
	void Output_field_k_inloop(IncVF& W);
	void Output_field_k_inloop(IncVF& W, IncSF& T);
	void Output_field_k_inloop(IncSF& T, string Pr_switch);
	void Output_field_k_inloop(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_pressure_spectrum_inloop();
	
	
	void Output_global();
	void Output_global(IncSF& T);
	void Output_global(IncVF& W);
	void Output_global(IncVF& W, IncSF& T);
	void Output_global(IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling);
	void Output_global(IncVF& W, IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling);
  
	void Output_field();
	void Output_field(IncSF& T);
	void Output_field(IncVF& W);
	void Output_field(IncVF& W, IncSF& T);
	void Output_field(IncSF& T, string Pr_switch);
	void Output_field(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_field_hdf5();

	void Output_field_frequent();
	void Output_field_frequent(IncSF& T);
	void Output_field_frequent(IncVF& W);
	void Output_field_frequent(IncVF& W, IncSF& T);
	void Output_field_frequent(IncSF& T, string Pr_switch);
	void Output_field_frequent(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_realfield();
	void Output_realfield(IncSF& T);
	void Output_realfield(IncVF& W);
	void Output_realfield(IncVF& W, IncSF& T);
	void Output_realfield(IncSF& T, string Pr_switch);
	void Output_realfield(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_realfield_hdf5();

	void Output_shell_spectrum();
	void Output_shell_spectrum(IncSF& T);
	void Output_shell_spectrum(IncVF& W);
	void Output_shell_spectrum(IncVF& W, IncSF &T);
	void Output_shell_spectrum(IncSF &T, string Pr_switch);
	void Output_shell_spectrum(IncVF& W, IncSF &T, string Pr_switch);
	
	void Output_ring_spectrum();
	void Output_ring_spectrum(IncSF& T);
	void Output_ring_spectrum(IncVF& W);
	void Output_ring_spectrum(IncVF& W, IncSF &T);
	void Output_ring_spectrum(IncSF &T, string Pr_switch);
	void Output_ring_spectrum(IncVF& W, IncSF &T, string Pr_switch);

	void Output_cylinder_ring_spectrum();
	void Output_cylinder_ring_spectrum(IncSF& T);
	void Output_cylinder_ring_spectrum(IncVF& W);
	void Output_cylinder_ring_spectrum(IncVF& W, IncSF &T);
	void Output_cylinder_ring_spectrum(IncSF &T, string Pr_switch);
	void Output_cylinder_ring_spectrum(IncVF& W, IncSF &T, string Pr_switch);
	
	void Output_pressure_spectrum();
	
	
	void Output_flux();
	void Output_flux(IncSF& T);
	void Output_flux(IncVF& W);
	void Output_flux(IncVF& W, IncSF& T);	
	void Output_flux(IncSF& T, string Pr_switch);
	void Output_flux(IncVF& W, IncSF &T, string Pr_switch);
	
	
	void Output_shell_to_shell();
	void Output_shell_to_shell(IncSF& T);
	void Output_shell_to_shell(IncVF& W);
	void Output_shell_to_shell(IncVF& W, IncSF& T);
	void Output_shell_to_shell(IncSF& T, string Pr_switch);
	void Output_shell_to_shell(IncVF& W, IncSF& T, string Pr_switch);

	void Output_field_reduced();
	void Output_field_reduced(IncSF& T);
	void Output_field_reduced(IncVF& W);
	void Output_field_reduced(IncVF& W, IncSF& T);
	void Output_field_reduced(IncSF& T, string Pr_switch);
	void Output_field_reduced(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_field_k();
	void Output_field_k(IncSF& T);
	void Output_field_k(IncVF& W);
    void Output_field_k(IncVF& W, IncSF& T);
	void Output_field_k(IncSF& T, string Pr_switch);
	void Output_field_k(IncVF& W, IncSF& T, string Pr_switch);

	void Output_field_r();
	void Output_field_r(IncSF& T);
	void Output_field_r(IncVF& W);
    void Output_field_r(IncVF& W, IncSF& T);
	void Output_field_r(IncSF& T, string Pr_switch);
	void Output_field_r(IncVF& W, IncSF& T, string Pr_switch);
	
	
	void Output_cout();
	void Output_cout(IncSF& T);
	void Output_cout(IncVF& W);
	void Output_cout(IncVF& W, IncSF& T);
	
	
	void Output_ring_to_ring();
	void Output_ring_to_ring(IncSF& T);
	void Output_ring_to_ring(IncVF& W);
	void Output_ring_to_ring(IncVF& W, IncSF &T);
	void Output_ring_to_ring(IncSF& T, string Pr_switch);
	void Output_ring_to_ring(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_cylinder_ring_to_ring();
	void Output_cylinder_ring_to_ring(IncSF& T);
	void Output_cylinder_ring_to_ring(IncVF& W);
	void Output_cylinder_ring_to_ring(IncVF& W, IncSF &T);
	void Output_cylinder_ring_to_ring(IncSF& T, string Pr_switch);
	void Output_cylinder_ring_to_ring(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_structure_fn();
	void Output_structure_fn(IncSF& T);
	void Output_structure_fn(IncVF& W);
	void Output_structure_fn(IncVF& W, IncSF &T);
	void Output_structure_fn(IncSF& T, string Pr_switch);
	void Output_structure_fn(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_planar_structure_fn();
	void Output_planar_structure_fn(IncSF& T);
	void Output_planar_structure_fn(IncVF& W);
	void Output_planar_structure_fn(IncVF& W, IncSF &T);
	void Output_planar_structure_fn(IncSF& T, string Pr_switch);
	void Output_planar_structure_fn(IncVF& W, IncSF& T, string Pr_switch);
	
	void Output_moment();
	void Output_moment(IncSF& T);
	void Output_moment(IncVF& W);
	void Output_moment(IncVF& W, IncSF &T);
	void Output_moment(IncSF& T, string Pr_switch);
	void Output_moment(IncVF& W, IncSF &T, string Pr_switch);
	
	void Output_Skpq();
	void Output_Skpq(IncSF& T);
	void Output_Skpq(IncVF& W);
	void Output_Skpq(IncVF& W, IncSF &T);
	void Output_Skpq(IncSF& T, string Pr_switch);
	void Output_Skpq(IncVF& W, IncSF &T, string Pr_switch);
	
	
};

#endif

//========================= Class declaration of IncFluid ends ============================== 
 

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

/*! \file  output_main.cc
 * 
 * @brief  Output_all_inloop & Output_field_k_inloop. 
 *
 *	In Output_all_inloop <BR>
 *	Output_global:	Total energy, dissipation etc. <BR>
 *	Output_cout:	Total energy, etc. on cout.	<BR>
 *	Output_field:	Field in k space. <BR?
 *	Output_field_reduced: Field in smaller box. <BR>
 *	Output_realfield: Field in r space <BR>
 *	Output_shell_spectrum: Shell spectrum <BR>
 *	Output_flux:	flux <BR>
 *	Output_shell_to_shell: shell-to-shell energy tr <BR>
 *	Output_field_frequent: frequent saving of field vars; Overwrites the field. <BR>
 *	Output_ring_spectrum: ring <BR>
 *	Output_ring_to_ring: ring-to-ring enregy transfer <BR>
 *	Output_cylinder_ring_spectrum: cylinderical ring spectrum <BR>
 *	Output_cylinder_ring_to_ring: Cylindrical ring-to-ring energy transfer <BR>
 *  Output_structire_fn: Structure function <BR>
 *	Output_planar_structure_fn: Planar structure function. <BR>
 *
 *	In Output_field_k_inloop: <BR>
 *	Output_field_k:	Outputs \f$ V(k) \f$ at specified k's.
 *	Output_field_r: Outputs \f$ V(r) \f$ at specified r's.
 *	Output_Skpq: Outputs S(k|p|q) for specified triads.
 *
 *	Output_pressure_spectrum_inloop: Outputs pressure spectrum. <BR>
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 


#include "../IncFluid.h"   


//*********************************************************************************************

/** @brief Outputs various things in main loop.
 *
 *	Outputs every save_next time. Saving takes place immediately after Tnow
 *	crosses save_next. <BR>
 *	save_next is updated to save_next + save_interval.
 *
 *	save_interval etc. are saved in Time class.
 */
/*
void IncFluid::Output_all_inloop_hdf5()
{
	if (Tnow >= Trealfield_save_next) 
	{
		Output_realfield_hdf5();
		Trealfield_save_next += Trealfield_save_interval;
		Close_output_files_hdf5("realField");
	}	

	if (Tnow >= Tfield_save_next)       
	{  
		Output_field_hdf5(); 
		Tfield_save_next += Tfield_save_interval;
		Close_output_files_hdf5("complexField");
	}	
	
	Close_output_files_hdf5("final");
}
*/
void IncFluid::Output_all_inloop()
{

	if (Tnow >= Tglobal_save_next) 
	{	 
		Output_global(); 
		Tglobal_save_next +=  Tglobal_save_interval; 
	}
			
	if (Tnow >= Tfield_save_next)       
	{  
		Output_field(); 
		Tfield_save_next += Tfield_save_interval;
	}	
	
	if (Tnow >= Tfield_frequent_save_next)       
	{  
		Output_field_frequent(); 
		Tfield_frequent_save_next += Tfield_frequent_save_interval;
	}	
	
	if (Tnow >= Tfield_reduced_save_next) 
	{	
		Output_field_reduced();
		Tfield_reduced_save_next += Tfield_reduced_save_interval;
	}	
	
	if (Tnow >= Trealfield_save_next) 
	{
		Output_realfield();
		Trealfield_save_next += Trealfield_save_interval;
	}	
			
	if (Tnow >= Tspectrum_save_next)
	{
		Output_shell_spectrum();	
		Tspectrum_save_next += Tspectrum_save_interval;
	}		
					
	if (Tnow >= Tflux_save_next)
	{
		Output_flux();			// 0 -> real_imag off
		Tflux_save_next += Tflux_save_interval;
	}	
			
	if (Tnow >= Tshell_to_shell_save_next)
	{
		Output_shell_to_shell();
		Tshell_to_shell_save_next += Tshell_to_shell_save_interval;
	}	
	

	if ((Tnow >= Tring_spectrum_save_next) && (CV_anisotropic_ring_switch == 1))
	{  
		Output_ring_spectrum(); 
		Tring_spectrum_save_next += Tring_spectrum_save_interval;
	}	
	
	if ((Tnow >= Tring_to_ring_save_next) && (CV_anisotropic_ring_switch == 1)) 
	{  
		Output_ring_to_ring(); 
		Tring_to_ring_save_next += Tring_to_ring_save_interval;
	}	
	
	if ((Tnow >= Tcylinder_ring_spectrum_save_next) && (CV_anisotropic_cylinder_switch == 1))
	{  
		Output_cylinder_ring_spectrum(); 
		Tcylinder_ring_spectrum_save_next += Tcylinder_ring_spectrum_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_to_ring_save_next) && (CV_anisotropic_cylinder_switch == 1)) 
	{  
		Output_cylinder_ring_to_ring(); 
		Tcylinder_ring_to_ring_save_next += Tcylinder_ring_to_ring_save_interval;
	}
	
	if (Tnow >= Tstructure_fn_save_next)
	{
		if (CV_structure_fn_switch == 1)
			Output_structure_fn();			
		
		if (CV_planar_structure_fn_switch == 1)
			Output_planar_structure_fn();
			
		Tstructure_fn_save_next += Tstructure_fn_save_interval;
	}
	
	if (Tnow >= Tmoment_save_next)
	{
		Output_moment();			
		Tmoment_save_next += Tmoment_save_interval;
	}
	
	if (Tnow >= Tcout_save_next)   
	{  
		Output_cout(); 
		Tcout_save_next += Tcout_save_interval; 
	}
}


//*********************************************************************************************
// with scalar
//
void IncFluid::Output_all_inloop(IncSF& T)
{
	
	if (Tnow >= Tglobal_save_next) 
	{	 
		Output_global(T); 
		Tglobal_save_next +=  Tglobal_save_interval; 
	}
	
			
	if (Tnow >= Tfield_save_next)       
	{  
		Output_field(T); 
		Tfield_save_next += Tfield_save_interval;
	}	
	
	if (Tnow >= Tfield_frequent_save_next)       
	{  
		Output_field_frequent(T); 
		Tfield_frequent_save_next += Tfield_frequent_save_interval;
	}
	
	
	if (Tnow >= Tfield_reduced_save_next) 
	{	
		Output_field_reduced(T);
		Tfield_reduced_save_next += Tfield_reduced_save_interval;
	}	
	
	if (Tnow >= Trealfield_save_next) 
	{
		Output_realfield(T);
		Trealfield_save_next += Trealfield_save_interval;
	}	
			
	if (Tnow >= Tspectrum_save_next)
	{
		Output_shell_spectrum(T);	
		Tspectrum_save_next += Tspectrum_save_interval;
	}		
			
	if (Tnow >= Tflux_save_next)
	{
		Output_flux(T);			// 0 -> real_imag off
		Tflux_save_next += Tflux_save_interval;
	}	
			
	if (Tnow >= Tshell_to_shell_save_next)
	{
		Output_shell_to_shell(T);
		Tshell_to_shell_save_next += Tshell_to_shell_save_interval;
	}	
		

	if ((Tnow >= Tring_spectrum_save_next)  && (CV_anisotropic_ring_switch == 1))
	{  
		Output_ring_spectrum(T); 
		Tring_spectrum_save_next += Tring_spectrum_save_interval;
	}	
	
	if ((Tnow >= Tring_to_ring_save_next) && (CV_anisotropic_ring_switch == 1))      
	{  
		Output_ring_to_ring(T); 
		Tring_to_ring_save_next += Tring_to_ring_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_spectrum_save_next) && (CV_anisotropic_cylinder_switch == 1))
	{  
		Output_cylinder_ring_spectrum(T); 
		Tcylinder_ring_spectrum_save_next += Tcylinder_ring_spectrum_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_to_ring_save_next) && (CV_anisotropic_cylinder_switch == 1)) 
	{  
		Output_cylinder_ring_to_ring(T); 
		Tcylinder_ring_to_ring_save_next += Tcylinder_ring_to_ring_save_interval;
	}
	
	if (Tnow >= Tstructure_fn_save_next)
	{
		if (CV_structure_fn_switch == 1)
			Output_structure_fn(T);			
		
		if (CV_planar_structure_fn_switch == 1)
			Output_planar_structure_fn(T);
			
		Tstructure_fn_save_next += Tstructure_fn_save_interval;
	}
	
	if (Tnow >= Tmoment_save_next)
	{
		Output_moment(T);			
		Tmoment_save_next += Tmoment_save_interval;
	}
	
	if (Tnow >= Tcout_save_next)   
	{  
		Output_cout(T); 
		Tcout_save_next += Tcout_save_interval; 
	}
}


//*********************************************************************************************
// with vector
//
void IncFluid::Output_all_inloop(IncVF& W)
{

	if (Tnow >= Tglobal_save_next) 
	{	 
		Output_global(W); 
		Tglobal_save_next +=  Tglobal_save_interval; 
	}
	
		
	if (Tnow >= Tfield_save_next)       
	{  
		Output_field(W); 
		Tfield_save_next += Tfield_save_interval;
	}	
	
	if (Tnow >= Tfield_frequent_save_next)       
	{  
		Output_field_frequent(W); 
		Tfield_frequent_save_next += Tfield_frequent_save_interval;
	}
	
	
	if (Tnow >= Tfield_reduced_save_next) 
	{	
		Output_field_reduced(W);
		Tfield_reduced_save_next += Tfield_reduced_save_interval;
	}	
	
	if (Tnow >= Trealfield_save_next) 
	{
		Output_realfield(W);
		Trealfield_save_next += Trealfield_save_interval;
	}	
	
			
	if (Tnow >= Tspectrum_save_next)
	{
		Output_shell_spectrum(W);	
		Tspectrum_save_next += Tspectrum_save_interval;
	}		
	
				
	if (Tnow >= Tflux_save_next)
	{
		Output_flux(W);			// 0 -> real_imag off
		Tflux_save_next += Tflux_save_interval;
	}	
			
			
		
	if (Tnow >= Tshell_to_shell_save_next)
	{
		Output_shell_to_shell(W);
		Tshell_to_shell_save_next += Tshell_to_shell_save_interval;
	}	
	

	if ((Tnow >= Tring_spectrum_save_next) && (CV_anisotropic_ring_switch == 1))       
	{  
		Output_ring_spectrum(W); 
		Tring_spectrum_save_next += Tring_spectrum_save_interval;
	}	
	
	if ((Tnow >= Tring_to_ring_save_next) && (CV_anisotropic_ring_switch == 1))      
	{  
		Output_ring_to_ring(W); 
		Tring_to_ring_save_next += Tring_to_ring_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_spectrum_save_next) && (CV_anisotropic_cylinder_switch == 1))
	{  
		Output_cylinder_ring_spectrum(W); 
		Tcylinder_ring_spectrum_save_next += Tcylinder_ring_spectrum_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_to_ring_save_next) && (CV_anisotropic_cylinder_switch == 1)) 
	{  
		Output_cylinder_ring_to_ring(W); 
		Tcylinder_ring_to_ring_save_next += Tcylinder_ring_to_ring_save_interval;
	}
	
		
	if (Tnow >= Tstructure_fn_save_next)
	{
		if (CV_structure_fn_switch == 1)
			Output_structure_fn(W);			
		
		if (CV_planar_structure_fn_switch == 1)
			Output_planar_structure_fn(W);
			
		Tstructure_fn_save_next += Tstructure_fn_save_interval;
	}
	
		
	if (Tnow >= Tmoment_save_next)
	{
		Output_moment(W);			
		Tmoment_save_next += Tmoment_save_interval;
	}
	
		
	if (Tnow >= Tcout_save_next)   
	{  
		Output_cout(W); 
		Tcout_save_next += Tcout_save_interval; 
	}
		
}




//*********************************************************************************************
// W and T
//
void IncFluid::Output_all_inloop(IncVF& W, IncSF& T)
{

	if (Tnow >= Tglobal_save_next) 
	{	 
		Output_global(W, T); 
		Tglobal_save_next +=  Tglobal_save_interval; 
	}
			
	if (Tnow >= Tfield_save_next)       
	{  
		Output_field(W, T); 
		Tfield_save_next += Tfield_save_interval;
	}	
	
	if (Tnow >= Tfield_frequent_save_next)       
	{  
		Output_field_frequent(W, T); 
		Tfield_frequent_save_next += Tfield_frequent_save_interval;
	}
	
	
	if (Tnow >= Tfield_reduced_save_next) 
	{	
		Output_field_reduced(W, T);
		Tfield_reduced_save_next += Tfield_reduced_save_interval;
	}	
	
	if (Tnow >= Trealfield_save_next) 
	{
		Output_realfield(W, T);
		Trealfield_save_next += Trealfield_save_interval;
	}	
			
	if (Tnow >= Tspectrum_save_next)
	{
		Output_shell_spectrum(W, T);	
		Tspectrum_save_next += Tspectrum_save_interval;
	}		
			
	if (Tnow >= Tflux_save_next)
	{
		Output_flux(W, T);			// 0 -> real_imag off
		Tflux_save_next += Tflux_save_interval;
	}	
			
	if (Tnow >= Tshell_to_shell_save_next)
	{
		Output_shell_to_shell(W, T);
		Tshell_to_shell_save_next += Tshell_to_shell_save_interval;
	}	

	if ((Tnow >= Tring_spectrum_save_next)  && (CV_anisotropic_ring_switch == 1))     
	{  
		Output_ring_spectrum(W, T); 
		Tring_spectrum_save_next += Tring_spectrum_save_interval;
	}	
	
	if ((Tnow >= Tring_to_ring_save_next) && (CV_anisotropic_ring_switch == 1))
	{  
		Output_ring_to_ring(W, T); 
		Tring_to_ring_save_next += Tring_to_ring_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_spectrum_save_next) && (CV_anisotropic_cylinder_switch == 1))
	{  
		Output_cylinder_ring_spectrum(W, T); 
		Tcylinder_ring_spectrum_save_next += Tcylinder_ring_spectrum_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_to_ring_save_next) && (CV_anisotropic_cylinder_switch == 1)) 
	{  
		Output_cylinder_ring_to_ring(W, T); 
		Tcylinder_ring_to_ring_save_next += Tcylinder_ring_to_ring_save_interval;
	}
	
	if (Tnow >= Tstructure_fn_save_next)
	{
		if (CV_structure_fn_switch == 1)
			Output_structure_fn(W, T);			
		
		if (CV_planar_structure_fn_switch == 1)
			Output_planar_structure_fn(W, T);
			
		Tstructure_fn_save_next += Tstructure_fn_save_interval;
	}
	
	if (Tnow >= Tmoment_save_next)
	{
		Output_moment(W, T);			
		Tmoment_save_next += Tmoment_save_interval;
	}
	
	if (Tnow >= Tcout_save_next)   
	{  
		Output_cout(W, T); 
		Tcout_save_next += Tcout_save_interval; 
	}
}



//*********************************************************************************************
// for convection
//
void IncFluid::Output_all_inloop(IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling)
{
	
	if (Tnow >= Tglobal_save_next) 
	{	 
		Output_global(T, Ra, Pr, Pr_switch, RB_Uscaling); 
		Tglobal_save_next +=  Tglobal_save_interval; 
	}
			
	if (Tnow >= Tfield_save_next)       
	{  
		Output_field(T, Pr_switch); 
		Tfield_save_next += Tfield_save_interval;
	}	
	
	if (Tnow >= Tfield_frequent_save_next)       
	{  
		Output_field_frequent(T, Pr_switch); 
		Tfield_frequent_save_next += Tfield_frequent_save_interval;
	}
	
	if (Tnow >= Tfield_reduced_save_next) 
	{	
		Output_field_reduced(T, Pr_switch);
		Tfield_reduced_save_next += Tfield_reduced_save_interval;
	}	
	
	if (Tnow >= Trealfield_save_next) 
	{
		Output_realfield(T, Pr_switch);
		Trealfield_save_next += Trealfield_save_interval;
	}	
			
	if (Tnow >= Tspectrum_save_next)
	{
		Output_shell_spectrum(T, Pr_switch);	
		Tspectrum_save_next += Tspectrum_save_interval;
	}		
			
	if (Tnow >= Tflux_save_next)
	{
		Output_flux(T, Pr_switch);			// 0 -> real_imag off
		Tflux_save_next += Tflux_save_interval;
	}	
			
	if (Tnow >= Tshell_to_shell_save_next)
	{
		Output_shell_to_shell(T, Pr_switch);
		Tshell_to_shell_save_next += Tshell_to_shell_save_interval;
	}	
	
	if ((Tnow >= Tring_spectrum_save_next) && (CV_anisotropic_ring_switch == 1))
	{  
		Output_ring_spectrum(T, Pr_switch); 
		Tring_spectrum_save_next += Tring_spectrum_save_interval;
	}	
	
	if ((Tnow >= Tring_to_ring_save_next) && (CV_anisotropic_ring_switch == 1))
	{  
		Output_ring_to_ring(T, Pr_switch); 
		Tring_to_ring_save_next += Tring_to_ring_save_interval;
	}	
	
	if ((Tnow >= Tcylinder_ring_spectrum_save_next) && (CV_anisotropic_cylinder_switch == 1))
	{  
		Output_cylinder_ring_spectrum(T, Pr_switch); 
		Tcylinder_ring_spectrum_save_next += Tcylinder_ring_spectrum_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_to_ring_save_next) && (CV_anisotropic_cylinder_switch == 1)) 
	{  
		Output_cylinder_ring_to_ring(T, Pr_switch); 
		Tcylinder_ring_to_ring_save_next += Tcylinder_ring_to_ring_save_interval;
	}
	
	if (Tnow >= Tstructure_fn_save_next)
	{
		if (CV_structure_fn_switch == 1)
			Output_structure_fn(T, Pr_switch);			
		
		if (CV_planar_structure_fn_switch == 1)
			Output_planar_structure_fn(T, Pr_switch);
			
		Tstructure_fn_save_next += Tstructure_fn_save_interval;
	}	
	
	if (Tnow >= Tmoment_save_next)
	{
		Output_moment(T, Pr_switch);			
		Tmoment_save_next += Tmoment_save_interval;
	}
	
	if (Tnow >= Tcout_save_next)   
	{  
		Output_cout(T); 
		Tcout_save_next += Tcout_save_interval; 
	}
					
}


//*********************************************************************************************
// Magnetoconvection
//
void IncFluid::Output_all_inloop(IncVF& W, IncSF& T,  DP Ra, DP Pr, string Pr_switch, 
									string RB_Uscaling)
{
		if (Tnow >= Tglobal_save_next) 
	{	 
		Output_global(W, T, Ra, Pr, Pr_switch, RB_Uscaling);
		Tglobal_save_next +=  Tglobal_save_interval; 
	}
			
	if (Tnow >= Tfield_save_next)       
	{  
		Output_field(W, T, Pr_switch); 
		Tfield_save_next += Tfield_save_interval;
	}	
	
	if (Tnow >= Tfield_frequent_save_next)       
	{  
		Output_field_frequent(W, T, Pr_switch); 
		Tfield_frequent_save_next += Tfield_frequent_save_interval;
	}
	
	if (Tnow >= Tfield_reduced_save_next) 
	{	
		Output_field_reduced(W, T, Pr_switch);
		Tfield_reduced_save_next += Tfield_reduced_save_interval;
	}	
	
	if (Tnow >= Trealfield_save_next) 
	{
		Output_realfield(W, T, Pr_switch);
		Trealfield_save_next += Trealfield_save_interval;
	}	
			
	if (Tnow >= Tspectrum_save_next)
	{
		Output_shell_spectrum(W, T, Pr_switch);	
		Tspectrum_save_next += Tspectrum_save_interval;
	}		
			
	if (Tnow >= Tflux_save_next)
	{
		Output_flux(W, T, Pr_switch);			// 0 -> real_imag off
		Tflux_save_next += Tflux_save_interval;
	}	
			
	if (Tnow >= Tshell_to_shell_save_next)
	{
		Output_shell_to_shell(W, T, Pr_switch);
		Tshell_to_shell_save_next += Tshell_to_shell_save_interval;
	}	
	
	if ((Tnow >= Tring_spectrum_save_next)  && (CV_anisotropic_ring_switch == 1))
	{  
		Output_ring_spectrum(W, T, Pr_switch); 
		Tring_spectrum_save_next += Tring_spectrum_save_interval;
	}	
	
	if ((Tnow >= Tring_to_ring_save_next) && (CV_anisotropic_ring_switch == 1))
	{  
		Output_ring_to_ring(W, T, Pr_switch); 
		Tring_to_ring_save_next += Tring_to_ring_save_interval;
	}	
	
	if ((Tnow >= Tcylinder_ring_spectrum_save_next) && (CV_anisotropic_cylinder_switch == 1))
	{  
		Output_cylinder_ring_spectrum(W, T, Pr_switch); 
		Tcylinder_ring_spectrum_save_next += Tcylinder_ring_spectrum_save_interval;
	}
	
	if ((Tnow >= Tcylinder_ring_to_ring_save_next) && (CV_anisotropic_cylinder_switch == 1)) 
	{  
		Output_cylinder_ring_to_ring(W, T, Pr_switch); 
		Tcylinder_ring_to_ring_save_next += Tcylinder_ring_to_ring_save_interval;
	}
	
	if (Tnow >= Tstructure_fn_save_next)
	{
		if (CV_structure_fn_switch == 1)
			Output_structure_fn(W, T, Pr_switch);			
		
		if (CV_planar_structure_fn_switch == 1)
			Output_planar_structure_fn(W, T, Pr_switch);
			
		Tstructure_fn_save_next += Tstructure_fn_save_interval;
	}	
	
	if (Tnow >= Tmoment_save_next)
	{
		Output_moment(W, T, Pr_switch);			
		Tmoment_save_next += Tmoment_save_interval;
	}	
	
	
	if (Tnow >= Tcout_save_next)   
	{  
		Output_cout(W, T); 
		Tcout_save_next += Tcout_save_interval; 
	}																						
}


//*********************************************************************************************

// Outputs V(k) at specified k's
void IncFluid::Output_field_k_inloop()
{
	if (Tnow >= Tfield_k_save_next)
	{
		Output_field_k();
		
		Tfield_k_save_next += Tfield_k_save_interval;
	}
	
	
	if (Tnow >= Tfield_r_save_next)
	{
		if (output_field_r_switch == 1)
			Output_field_r();
			
		Tfield_r_save_next += Tfield_r_save_interval;
	}
	
	
	if (Tnow >= TSkpq_save_next)
	{
		if (skpq_switch	== 1)
			Output_Skpq();
			
		TSkpq_save_next += TSkpq_save_interval;
	}
}

void IncFluid::Output_field_k_inloop(IncSF& T)
{
	if (Tnow >= Tfield_k_save_next)
	{
		Output_field_k(T);
		
		Tfield_k_save_next += Tfield_k_save_interval;
	}
	
	if (Tnow >= Tfield_r_save_next)
	{
		if (output_field_r_switch == 1) 
			Output_field_r(T);
			
		Tfield_r_save_next += Tfield_r_save_interval;
	}
	
	
	if (Tnow >= TSkpq_save_next)
	{
		if (skpq_switch	== 1)
			Output_Skpq(T);
			
		TSkpq_save_next += TSkpq_save_interval;
	}
	
}


void IncFluid::Output_field_k_inloop(IncVF& W)
{
	if (Tnow >= Tfield_k_save_next)
	{
		Output_field_k(W);
		
		Tfield_k_save_next += Tfield_k_save_interval;
	}
	
	
	if (Tnow >= Tfield_r_save_next)
	{
		if (output_field_r_switch == 1)
			Output_field_r(W);
			
		Tfield_r_save_next += Tfield_r_save_interval;
	}
	
	
	if (Tnow >= TSkpq_save_next)
	{
		if (skpq_switch	== 1)
			Output_Skpq(W);
			
		TSkpq_save_next += TSkpq_save_interval;
	}
}

void IncFluid::Output_field_k_inloop(IncVF& W, IncSF& T)
{
	if (Tnow >= Tfield_k_save_next)
	{
		Output_field_k(W, T);
		
		Tfield_k_save_next += Tfield_k_save_interval;
	}
	
	
	if (Tnow >= Tfield_r_save_next)
	{
		if (output_field_r_switch == 1)
			Output_field_r(W, T);
			
		Tfield_r_save_next += Tfield_r_save_interval;
	}
	
	
	if (Tnow >= TSkpq_save_next)
	{
		if (skpq_switch	== 1)
			Output_Skpq(W, T);
			
		TSkpq_save_next += TSkpq_save_interval;
	}
}



void IncFluid::Output_field_k_inloop(IncSF& T, string Pr_switch)
{

	if (Tnow >= Tfield_k_save_next)
	{
		Output_field_k(T, Pr_switch);
		
		Tfield_k_save_next += Tfield_k_save_interval;
	}
	
	if (Tnow >= Tfield_r_save_next)
	{
		if (output_field_r_switch == 1)
			Output_field_r(T, Pr_switch);
			
		Tfield_r_save_next += Tfield_r_save_interval;
	}
	
	
	if (Tnow >= TSkpq_save_next)
	{
		if (skpq_switch	== 1)
			Output_Skpq(T, Pr_switch);
			
		TSkpq_save_next += TSkpq_save_interval;
	}
}



void IncFluid::Output_field_k_inloop(IncVF& W, IncSF& T, string Pr_switch)
{
	if (Tnow >= Tfield_k_save_next)
	{
		Output_field_k(W, T, Pr_switch);
		
		Tfield_k_save_next += Tfield_k_save_interval;
	}
	
	
	if (Tnow >= Tfield_r_save_next)
	{
		if (output_field_r_switch == 1)
			Output_field_r(W, T, Pr_switch);
			
		Tfield_r_save_next += Tfield_r_save_interval;
	}
	
	
	if (Tnow >= TSkpq_save_next)
	{
		if (skpq_switch	== 1)
			Output_Skpq(W, T, Pr_switch);
			
		TSkpq_save_next += TSkpq_save_interval;
	}
}



//*********************************************************************************************

// Outputs pressure spectrum p(k)
void IncFluid::Output_pressure_spectrum_inloop()
{
	if (Tnow >= Tspectrum_pressure_save_next)
	{
		Output_pressure_spectrum();	
		
		Tspectrum_pressure_save_next += Tspectrum_pressure_save_interval;
	}		
}




//*********************************************************************************************

// Outputs Tnow and total energy to cout.
void IncFluid::Output_cout()
{	
	CV_Compute_totalenergy_diss();
	
	if (my_id == master_id)
		cout  << Tnow << " "  << CV_total_energy << endl;
}


void IncFluid::Output_cout(IncSF& T)
{
	CV_Compute_totalenergy_diss(); 
	T.CS_Compute_totalenergy_diss();
	
	if (my_id == master_id)
		cout  << Tnow << " "  << CV_total_energy << " " << T.CS_total_energy  << endl;
}

void IncFluid::Output_cout(IncVF& W)
{
	CV_Compute_totalenergy_diss(); 
	W.CV_Compute_totalenergy_diss();
		
	if (my_id == master_id)
		cout  << Tnow << " "  << CV_total_energy << " " << W.CV_total_energy << endl;
}

void IncFluid::Output_cout(IncVF& W, IncSF& T)		// for RB convection
{
	CV_Compute_totalenergy_diss(); 
	W.CV_Compute_totalenergy_diss();
	T.CS_Compute_totalenergy_diss();
	
	if (my_id == master_id)
		cout  << Tnow << " "  << CV_total_energy << " " << W.CV_total_energy << " " 
		 << T.CS_total_energy << endl;
}


//********************************  End of output_main.cc ************************************




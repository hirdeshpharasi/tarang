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

/*! \file  IncFluid_Time.cc
 * 
 * @brief  Class constructor of IncFluid_Time.   
 *
 * @sa   IncFluid_Time.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 
 
#include "IncFluid_Time.h"  


//********************************************************************************************* 

Time::Time(Array<DP,1> time_para, Array<DP,1> time_save_interval)   // time initialization
{
	
	Tinit				= time_para(1);
	Tfinal				= time_para(2);
	Tdt_fixed			= time_para(3);
	Tdiagnostics_init	= time_para(4);
	Tstructure_fn_start = time_para(5);
	
	Tnow = Tinit;
  
	Tglobal_save_interval			= time_save_interval(1);
	Tfield_save_interval			= time_save_interval(2);
	Tfield_frequent_save_interval	= time_save_interval(3);
	Tfield_reduced_save_interval	= time_save_interval(4);
	Trealfield_save_interval		= time_save_interval(5);
	Tfield_k_save_interval			= time_save_interval(6);
	Tfield_r_save_interval			= time_save_interval(7);
	Tspectrum_save_interval			= time_save_interval(8);
	Tspectrum_pressure_save_interval= time_save_interval(9);	
	Tflux_save_interval				= time_save_interval(10);
	Tshell_to_shell_save_interval	= time_save_interval(11);
	Tring_spectrum_save_interval	= time_save_interval(12);
	Tring_to_ring_save_interval		= time_save_interval(13);
	Tcylinder_ring_spectrum_save_interval =	time_save_interval(14);
	Tcylinder_ring_to_ring_save_interval  =	time_save_interval(15);
	Tstructure_fn_save_interval		= time_save_interval(16);
	Tmoment_save_interval			= time_save_interval(17);
	TSkpq_save_interval				= time_save_interval(18);		
	Tcout_save_interval				= time_save_interval(19);  
	
	Tglobal_save_next				= Tinit;	
	Tfield_save_next				= Tinit;
	Tfield_frequent_save_next		= Tinit + Tfield_frequent_save_interval; 
	Tfield_reduced_save_next		= Tinit;
	Trealfield_save_next			= Tinit;
	Tfield_k_save_next				= Tinit;
	Tfield_r_save_next				= Tdiagnostics_init;
	Tspectrum_save_next				= Tinit;
	Tspectrum_pressure_save_next	= Tspectrum_save_next;
	Tflux_save_next					= Tdiagnostics_init;
	Tshell_to_shell_save_next		= Tdiagnostics_init;	
	Tring_spectrum_save_next		= Tdiagnostics_init;
	Tring_to_ring_save_next			= Tdiagnostics_init;
	Tcylinder_ring_spectrum_save_next  = Tdiagnostics_init;
	Tcylinder_ring_to_ring_save_next   = Tdiagnostics_init;
	Tstructure_fn_save_next			= Tstructure_fn_start;
	Tmoment_save_next				= Tdiagnostics_init;
	TSkpq_save_next					= Tdiagnostics_init;
	Tcout_save_next					= Tinit;
}

//******************************* End of IncFluid_Time.cc  ************************************
  


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

/*! \file  IncFluid_Time.h
 * 
 * @brief  Class constructor of IncFluid_Time.   
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
              

#ifndef _H_IncFluid_Time
#define _H_IncFluid_Time

#include "../IncFlow/IncVF.h"

//*********************************************************************************************

class Time
{
	public:
	DP		Tinit; 
	DP		Tfinal;
	DP		Tdt_fixed;
	DP		Tdiagnostics_basic_init;						// Start basic diagnostics from
	DP		Tdiagnostics_advanced_init;						// Start advanced diagnostics from
	DP		Tstructure_fn_start;
	
	DP		Tdt;											// variable time including CFL
	DP		Tnow;
	
	DP		Tglobal_save_interval;	
	DP		Tfield_save_interval; 
	DP		Tfield_frequent_save_interval;
	DP		Tfield_reduced_save_interval;
	DP		Trealfield_save_interval;
	DP		Tfield_k_save_interval;				
	DP		Tfield_r_save_interval;	
	DP		Tspectrum_save_interval;
	DP		Tspectrum_pressure_save_interval;
	DP		Tflux_save_interval;
	DP		Tshell_to_shell_save_interval;		
	DP		Tring_spectrum_save_interval;
	DP		Tring_to_ring_save_interval;
	DP		Tcylinder_ring_spectrum_save_interval;
	DP		Tcylinder_ring_to_ring_save_interval;
	DP		Tstructure_fn_save_interval;
	DP		Tmoment_save_interval;
	DP		TSkpq_save_interval;						// Skpq saved here if skpq_switch =1
	DP		Tcout_save_interval;
	
	DP		Tglobal_save_next;	
	DP		Tfield_save_next; 
	DP		Tfield_frequent_save_next;
	DP		Tfield_reduced_save_next;
	DP		Trealfield_save_next;
	DP		Tfield_k_save_next;
	DP		Tfield_r_save_next;
	DP		Tspectrum_save_next;
	DP		Tspectrum_pressure_save_next;
	DP		Tflux_save_next;
	DP		Tshell_to_shell_save_next;	
	DP		Tring_spectrum_save_next;
	DP		Tring_to_ring_save_next;
	DP		Tcylinder_ring_spectrum_save_next;
	DP		Tcylinder_ring_to_ring_save_next;
	DP		Tstructure_fn_save_next;
	DP		Tmoment_save_next;
	DP		TSkpq_save_next;
	DP		Tcout_save_next;

	Time(Array<DP,1> time_para, Array<DP,1> time_save_interval);
};

#endif

//******************************* End of IncFluid_Time.h  *************************************
 

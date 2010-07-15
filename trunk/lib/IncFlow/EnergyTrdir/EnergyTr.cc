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
/*! \file EnergyTr.cc
 * 
 * @brief  Class declaration of EnergyTr 
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec 2008
 *
 @ @bug  No known bugs
 */
 
#include "../EnergyTr.h"


//*********************************************************************************************

EnergyTr::EnergyTr
(
	int *NN, 
	string string_switches[], 
	Array<int,1> switches, 
	DP *prog_kfactor, 
	Array<int,1> ET_parameters, 
	Array<DP,1> ET_shell_radii_sector_array
)						 
{
	int array_offset;
	
	
	ET_basis_type		= string_switches[1];
	ET_alias_switch		= string_switches[2];
	
	ET_anisotropic_ring_switch		= switches(4);
	ET_anisotropic_cylinder_switch	= switches(5);
	ET_Skpq_switch					= switches(8);
	
	// Extract shell-radii, sector etc. from ET_parameters()
	
	ET_real_imag_switch			= ET_parameters(1);
	no_spheres					= ET_parameters(2);
	no_shells					= ET_parameters(3);
	
	no_ring_shells				= ET_parameters(4);
	no_sectors_ring_tr			= ET_parameters(5);
	
	no_cylinder_shells			= ET_parameters(6);
	no_cylinder_kpll_slabs		= ET_parameters(7);
	
	int shell_input_scheme			= ET_parameters(8);
	int ring_shell_input_scheme		= ET_parameters(9);
	int sector_input_scheme			= ET_parameters(10);
	int cylinder_shell_input_scheme = ET_parameters(11);
	int cylinder_kpll_input_scheme  = ET_parameters(12);
	
	// Assign arrays

    for (int i=1; i<=3; i++) 
	{
      NET[i]=NN[i]; 
	  ET_kfactor[i] = prog_kfactor[i];
    } 
	
	 
	G1 = new Array<complx,3>(local_N1, NET[2], NET[3]/2+1);
	G2 = new Array<complx,3>(local_N1, NET[2], NET[3]/2+1);
	G3 = new Array<complx,3>(local_N1, NET[2], NET[3]/2+1);
	
	*G1 = 0.0; 
	*G2 = 0.0; 
	*G3 = 0.0; 


	
//************************  Initialize shell_radius and sphere_radius *************************
	
	DP max_radius_inside  
			= 1.0* Max_radius_inside(ET_basis_type, ET_alias_switch, NET, ET_kfactor);
	
	//
	// shell radii
	//
	
	if ( (no_shells < 6) || (max_radius_inside <= 16) ) 
	{   
		if (my_id == master_id)	
		{		
			cout << "WARNING: NO OF ENERGY TR SHELLS SHOULD BE GREATER THAN EQUAL TO SIX " 
				 << endl;
			cout <<	 "OR Maxpossible inner radius should be greater than equal to  16" << endl;
			cout << "Dont invoke flux and shell-to-shell routines. " << endl;
		}	
	}
	
	if (no_shells != no_spheres) 
	{
		cout << "ERROR: NO OF ENERGY TR SHELLS AND SPHERES ARE NOT EQUAL " << endl;
		cout << "EXITING RIGHT NOW.." << endl;
		exit(1);
	}
	
//************************  Initialize shell_radius and sphere_radius *************************	
	
	sphere_radius = new Array<DP,1>(no_shells+1);
	shell_radius = new Array<DP,1>(no_shells+1);
	
	temp_shell_tr = new Array<DP,1>(no_shells+1);
	temp_shell_tr_real = new Array<DP,1>(no_shells+1);
	temp_shell_tr_imag = new Array<DP,1>(no_shells+1);
	
	
	if (shell_input_scheme == 0)  // Fill using 2^(n/4) schele
	{
		(*shell_radius)(0) = (*sphere_radius)(0) = 0.0;  
		(*shell_radius)(1) = (*sphere_radius)(1) = 2.0; 
		(*shell_radius)(2) = (*sphere_radius)(2) = 4.0;  
		(*shell_radius)(3) = (*sphere_radius)(3) = 8.0;  
		
		(*shell_radius)(no_shells-2) = (*sphere_radius)(no_shells-2) 
									= max_radius_inside/2;
									
		(*shell_radius)(no_shells-1) = (*sphere_radius)(no_shells-1) 
									= max_radius_inside; 
																		
		(*shell_radius)(no_shells) = (*sphere_radius)(no_shells) = INF_RADIUS;
		
		// if no_shells==6, the above assignments will do (max_radius_inside > 16).
		
		if (no_shells > 6) 
		{		
			DP s = log2(max_radius_inside/16.0) / (no_shells - 5);
			
			for (int index = 4; index <= no_shells-3; index++) 
				(*shell_radius)(index) = (*sphere_radius)(index) = 8 * pow(2, s*(index-3));  	
		}		
	}
	
	else // use ET_shell_radii_sector_array() -- entries 0, r(1), r(2), .. r(n-1). Inf skipped.
	{ 
		for (int index = 0; index <= no_shells-1; index++) 
			(*shell_radius)(index) = (*sphere_radius)(index)  
								  = ET_shell_radii_sector_array(index); 	
		
		(*shell_radius)(no_shells) = (*sphere_radius)(no_shells) = INF_RADIUS;
											
	}	
	
	if (my_id == master_id)
		cout <<"ET: radii of the spheres are " << (*sphere_radius) << endl << endl;
	
	
	//************************   ring_shell radii *********************************************
	
	if (ET_anisotropic_ring_switch == 1)
	{
		if ( (no_ring_shells < 5) || (max_radius_inside <= 16) ) 
		{     
			if (my_id == master_id)	
			{
				cout << "WARNING: NO OF ENERGY TR RING_SHELLS SHOULD BE GREATER THAN"  
					 <<	"EQUAL TO FIVE"  << endl;
				cout <<	 "OR Maxpossible inner radius should be greater than equal to  16" 
					 << endl;
				cout << "Dont invoke flux and shell-to-shell routines. " << endl;
			}	
		}
		
		ring_shell_radius = new Array<DP,1>(no_ring_shells+1);
		
		
		temp_ring_tr = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
		temp_ring_tr_real = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
		temp_ring_tr_imag = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
		
		if (ring_shell_input_scheme == 0)  // Fill using 2^(n/4) scheme
		{
			(*ring_shell_radius)(0) = 0.0;  
			(*ring_shell_radius)(1) = 4.0; 
			(*ring_shell_radius)(2) = 8.0;  
			(*ring_shell_radius)(no_ring_shells-2) = max_radius_inside/2;
			(*ring_shell_radius)(no_ring_shells-1) = max_radius_inside; 
			(*ring_shell_radius)(no_ring_shells)   =  INF_RADIUS;
			// if no_ring_shells==5, the above assignments will do (max_radius_inside > 16).		
					
			if (no_ring_shells > 5) 
			{		
				DP s = log2(max_radius_inside/16.0) / (no_ring_shells - 4);
				
				for (int index = 3; index <= no_ring_shells-3; index++) 
					(*ring_shell_radius)(index) =  8*pow(2, s*(index-2));  	
			}		
		}
		
		else 	// use ET_shell_radii_sector_array() -- entries 0, r(1), r(2), .. r(n-1).
		{ 
		
			array_offset = no_shells;
			
			for (int index = 0; index<=no_ring_shells-1; index++) 
				(*ring_shell_radius)(index) = ET_shell_radii_sector_array(array_offset+index); 
														
			(*ring_shell_radius)(no_ring_shells) = INF_RADIUS;										
		}	
		
		if (my_id == master_id)	
			cout <<	"Radii of the rings:   "	<< (*ring_shell_radius) << endl << endl;
		
		//************************   Angles of the sectors ****************************************
		
		sector_angle_ring_tr = new Array<DP,1>(no_sectors_ring_tr+1);
		

		DP max_theta = Get_max_polar_angle(ET_basis_type);
		
		if (sector_input_scheme == 0)		// using linear spacing
		{
			DP dtheta = max_theta/no_sectors_ring_tr;
			
			for (int i=0; i<no_sectors_ring_tr; i++)
				(*sector_angle_ring_tr)(i) = i*dtheta;
					
			(*sector_angle_ring_tr)(no_sectors_ring_tr) = max_theta;
		}
		
		else if (sector_input_scheme == 1)		// same no of modes in each sector
		{
			;    // work on it
		}
		
		else if (sector_input_scheme == 2)	
		{
			array_offset  =  no_shells + no_ring_shells;
			
			for (int index=0; index < no_sectors_ring_tr; index++)
				(*sector_angle_ring_tr)(index) 
					= ET_shell_radii_sector_array(array_offset+index); 	
				
			(*sector_angle_ring_tr)(no_sectors_ring_tr) = max_theta;	
		}
																										

		if (my_id == master_id)	
			cout << "ring transfers: no_sectors & sector_angles:  " << no_sectors_ring_tr 
					<<  ":" << endl << (*sector_angle_ring_tr) << endl << endl;
								
	}   // END of if (ET_anisotropic_ring_switch == 1)
	
	
	//*****************************************************************************************
	
	// Now for cylinderical shells..
	//
	if (ET_anisotropic_cylinder_switch == 1)
	{		
		DP max_cylinder_radius_inside
			= 1.0* Anis_max_Krho_radius_inside(ET_basis_type, ET_alias_switch, 
													NET, ET_kfactor);
				
		
		if ( (no_cylinder_shells < 5) || (max_radius_inside <= 16) ) 
		{    
			if (my_id == master_id)
			{
				cout << "WARNING: NO OF ENERGY TR CYLINDER_SHELLS SHOULD BE GREATER " << 
						" THAN EQUAL TO FIVE" << endl;
				cout <<	 "OR Maxpossible inner radius should be greater than equal to  16" 
					 << endl;
				cout << "Dont invoke flux and shell-to-shell routines. " << endl;
			}	
		}

		
		cylinder_shell_radius = new Array<DP,1>(no_cylinder_shells+1);
		
		temp_cylinder_ring_tr = new Array<DP,2>(no_cylinder_shells+1, 
													no_cylinder_kpll_slabs+1);
													
		temp_cylinder_ring_tr_real = new Array<DP,2>(no_cylinder_shells+1, 
													 no_cylinder_kpll_slabs+1);
													 
		temp_cylinder_ring_tr_imag = new Array<DP,2>(no_cylinder_shells+1,  
													 no_cylinder_kpll_slabs+1);
		
		if (cylinder_shell_input_scheme == 0)  // Fill using 2^(n/4) scheme
		{
			(*cylinder_shell_radius)(0) = 0.0;  
			(*cylinder_shell_radius)(1) = 4.0; 
			(*cylinder_shell_radius)(2) = 8.0;  
			(*cylinder_shell_radius)(no_cylinder_shells-2) = max_cylinder_radius_inside/2;
			(*cylinder_shell_radius)(no_cylinder_shells-1) = max_cylinder_radius_inside; 
			(*cylinder_shell_radius)(no_cylinder_shells)   =  INF_RADIUS;
			// if no_cylindrical_shells==5, the above assignments will do 
			// (max_radius_inside > 16).
							
			if (no_ring_shells > 5) 
			{		
				DP s = log2(max_cylinder_radius_inside/16.0) / (no_cylinder_shells - 4);
				
				for (int index = 3; index <= no_cylinder_shells-3; index++) 
					(*cylinder_shell_radius)(index) =  8*pow(2, s*(index-2));  	
			}		
		}
		
		else 	// use ET_shell_radii_sector_array() -- entries 0, r(1), r(2), .. r(n-1).
		{ 	
			array_offset = no_shells+ no_ring_shells+ no_sectors_ring_tr;
			
			for (int index = 0; index<=no_cylinder_shells-1; index++) 
				(*cylinder_shell_radius)(index) 
						= ET_shell_radii_sector_array(array_offset+index); 
														
			(*cylinder_shell_radius)(no_cylinder_shells) = INF_RADIUS;										
		}	
		
		if (my_id == master_id)
			cout <<	"Radii of the cylinder rings:   "<< (*cylinder_shell_radius) 
				 << endl << endl;
		

		//   cylinder kpll assignment  
		//
		cylinder_kpll_array_tr = new Array<DP,1>(no_cylinder_kpll_slabs+1);
		
		DP	kkpll_min = 0.0;
		DP	kkpll_max = Anis_max_Kpll(ET_basis_type, ET_alias_switch, NET, ET_kfactor);
			
			
		if (cylinder_kpll_input_scheme == 0)
		{
			(*cylinder_kpll_array_tr)(0) = kkpll_min;
			
			DP dkkpll = (kkpll_max - kkpll_min)/ no_cylinder_kpll_slabs;
				
			for (int i=1; i<no_cylinder_kpll_slabs; i++)
				(*cylinder_kpll_array_tr)(i) = kkpll_min + i*dkkpll;
				
			(*cylinder_kpll_array_tr)(no_cylinder_kpll_slabs) = kkpll_max;
		}
		
		else 
		{
			array_offset = no_shells+ no_ring_shells+ no_sectors_ring_tr+ no_cylinder_shells;
			
			for (int index = 0; index<=no_cylinder_kpll_slabs-1; index++)
				(*cylinder_kpll_array_tr)(index) 
						= ET_shell_radii_sector_array(array_offset+index);
				
			(*cylinder_kpll_array_tr)(no_cylinder_kpll_slabs) = kkpll_max;	
		
		}
		
		if (my_id == master_id)
			cout << "cylindrical ring transfers: no_slabs & slab heights:  " 
				 << no_cylinder_kpll_slabs <<  ":" << endl 
				 << (*cylinder_kpll_array_tr) << endl << endl;

	}		// end of if (ET_anisotropic_cylinder_switch == 1)
	

	//************************   for S(k.p,q) calc		***************************************		
	
	if (ET_Skpq_switch == 1)
	{
		triad_array		= new Array<int,2>(MAX_NO_SKPQ_TRIADS+1, 6);
		*triad_array	= 0;
		
		Sself_array		= new Array<DP,1>(MAX_NO_SKPQ_TRIADS+1);
		S_to_VF_array	= new Array<DP,1>(MAX_NO_SKPQ_TRIADS+1);	
		S_SF_array		= new Array<DP,1>(MAX_NO_SKPQ_TRIADS+1);
		
		*Sself_array	= 0.0;
		*S_to_VF_array	= 0.0;
		*S_SF_array		= 0.0;
	}
	
	flux_self = new Array<DP,1>(no_shells+1);
	flux_VF_in_in = new Array<DP,1>(no_shells+1);
	flux_VF_in_out = new Array<DP,1>(no_shells+1);
	flux_VF_out_out = new Array<DP,1>(no_shells+1);
	flux_Elsasser = new Array<DP,1>(no_shells+1);
	flux_SF = new Array<DP,1>(no_shells+1);
	flux_hk = new Array<DP,1>(no_shells+1);
	
	*flux_self = 0.0;  
	*flux_VF_in_in = 0.0; 
	*flux_VF_in_out = 0.0; 
	*flux_Elsasser = 0.0;
	*flux_VF_out_out = 0.0; 
	*flux_SF = 0.0; 
	*flux_hk = 0.0;
	
	forceV_shell = new Array<DP,1>(no_shells+1);
	forceSF_shell = new Array<DP,1>(no_shells+1);
	*forceV_shell = 0.0; 
	*forceSF_shell = 0.0;
	
	
	shelltoshell_self = new Array<DP,2>(no_shells+1, no_shells+1);
	shelltoshell_VF = new Array<DP,2>(no_shells+1, no_shells+1);
	shelltoshell_Elsasser = new Array<DP,2>(no_shells+1, no_shells+1);
	shelltoshell_SF = new Array<DP,2>(no_shells+1, no_shells+1);
	shelltoshell_hk = new Array<DP,2>(no_shells+1, no_shells+1);
	
	*shelltoshell_self = 0.0;  
	*shelltoshell_VF = 0.0;
	*shelltoshell_SF = 0.0; 
	*shelltoshell_Elsasser = 0.0;
	*shelltoshell_hk = 0.0;
	
	energy_tr_shell_B0 =  new Array<DP,1>(no_shells+1);
	
	*energy_tr_shell_B0 = 0.0;
		
		
	if (ET_real_imag_switch == 1)
	{
		flux_self_real = new Array<DP,1>(no_shells+1);
		flux_VF_in_in_real = new Array<DP,1>(no_shells+1);
		flux_VF_in_out_real = new Array<DP,1>(no_shells+1);
		flux_VF_out_out_real = new Array<DP,1>(no_shells+1);
		flux_Elsasser_real = new Array<DP,1>(no_shells+1);
		flux_SF_real = new Array<DP,1>(no_shells+1);
		
		*flux_self_real = 0.0;  
		*flux_VF_in_in_real = 0.0; 
		*flux_VF_in_out_real = 0.0; 
		*flux_Elsasser_real = 0.0;
		*flux_VF_out_out_real = 0.0;
		*flux_SF_real = 0.0; 
		
		
		flux_self_imag = new Array<DP,1>(no_shells+1);
		flux_VF_in_in_imag = new Array<DP,1>(no_shells+1);
		flux_VF_in_out_imag = new Array<DP,1>(no_shells+1);
		flux_VF_out_out_imag = new Array<DP,1>(no_shells+1);
		flux_Elsasser_imag = new Array<DP,1>(no_shells+1);
		flux_SF_imag = new Array<DP,1>(no_shells+1);
		
		*flux_self_imag = 0.0;  
		*flux_VF_in_in_imag = 0.0; 
		*flux_VF_in_out_imag = 0.0; 
		*flux_Elsasser_imag = 0.0;
		*flux_VF_out_out_imag = 0.0; 
		*flux_SF_imag = 0.0; 
		
		forceV_shell_real = new Array<DP,1>(no_shells+1);
		forceSF_shell_real = new Array<DP,1>(no_shells+1);
		forceV_shell_imag = new Array<DP,1>(no_shells+1);
		forceSF_shell_imag = new Array<DP,1>(no_shells+1);
			
		*forceV_shell_real = 0.0;
		*forceSF_shell_real = 0.0;
		*forceV_shell_imag = 0.0; 
		*forceSF_shell_imag = 0.0;
		
		shelltoshell_self_real = new Array<DP,2>(no_shells+1, no_shells+1);
		shelltoshell_VF_real = new Array<DP,2>(no_shells+1, no_shells+1);
		shelltoshell_Elsasser_real = new Array<DP,2>(no_shells+1, no_shells+1);
		shelltoshell_SF_real = new Array<DP,2>(no_shells+1, no_shells+1);
		
		*shelltoshell_self_real = 0.0;  
		*shelltoshell_VF_real = 0.0;	
		*shelltoshell_SF_real = 0.0;
		*shelltoshell_Elsasser_real = 0.0;	
		
		shelltoshell_self_imag = new Array<DP,2>(no_shells+1, no_shells+1);
		shelltoshell_VF_imag = new Array<DP,2>(no_shells+1, no_shells+1);
		shelltoshell_Elsasser_imag = new Array<DP,2>(no_shells+1, no_shells+1);
		shelltoshell_SF_imag = new Array<DP,2>(no_shells+1, no_shells+1);
		
		*shelltoshell_self_imag = 0.0;  
		*shelltoshell_VF_imag = 0.0;	
		*shelltoshell_SF_imag = 0.0;
		*shelltoshell_Elsasser_imag = 0.0;
	}	
	
	// ring_to_ring
	
	
	if (ET_anisotropic_ring_switch == 1)
	{
		ring_to_ring_self		= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
									no_ring_shells+1, no_sectors_ring_tr+1);
									
		ring_to_ring_VF			= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
									no_ring_shells+1, no_sectors_ring_tr+1);
									
		ring_to_ring_Elsasser	= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
									no_ring_shells+1, no_sectors_ring_tr+1);
									
		ring_to_ring_SF			= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
									no_ring_shells+1, no_sectors_ring_tr+1);
		
		*ring_to_ring_self = 0.0;  
		*ring_to_ring_VF = 0.0;
		*ring_to_ring_SF = 0.0; 
		*ring_to_ring_Elsasser = 0.0;
		
		forceV_ring		= new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
		forceSF_ring	= new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
		
		*forceV_ring = 0.0;
		*forceSF_ring = 0.0;
		
		energy_tr_ring_B0  = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
		
		*energy_tr_ring_B0 = 0.0;
		
		
		if (ET_real_imag_switch == 1)
		{
			ring_to_ring_self_real	= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
										
			ring_to_ring_VF_real	= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
										
			ring_to_ring_Elsasser_real= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
										
			ring_to_ring_SF_real	= new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
			
			*ring_to_ring_self_real = 0.0;  
			*ring_to_ring_VF_real = 0.0;	
			*ring_to_ring_SF_real = 0.0;
			*ring_to_ring_Elsasser_real = 0.0;
			
			
			ring_to_ring_self_imag = new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
										
			ring_to_ring_VF_imag = new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
										
			ring_to_ring_Elsasser_imag = new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
										
			ring_to_ring_SF_imag = new Array<DP,4>(no_ring_shells+1, no_sectors_ring_tr+1, 
										no_ring_shells+1, no_sectors_ring_tr+1);
			
			*ring_to_ring_self_imag = 0.0;  
			*ring_to_ring_VF_imag = 0.0;	
			*ring_to_ring_SF_imag = 0.0;
			*ring_to_ring_Elsasser_imag = 0.0;
			
			
			forceV_ring_real = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
			forceSF_ring_real = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
			forceV_ring_imag = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
			forceSF_ring_imag = new Array<DP,2>(no_ring_shells+1, no_sectors_ring_tr+1);
			
			*forceV_ring_real = 0.0;
			*forceSF_ring_real = 0.0;
			*forceV_ring_imag = 0.0;
			*forceSF_ring_imag = 0.0;
		}
	}

	//
	// For cylinder
	//
	if (ET_anisotropic_cylinder_switch == 1)
	{	
		cylinder_ring_to_ring_self = new Array<DP,4>(no_cylinder_shells+1,
				no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);
		
		cylinder_ring_to_ring_VF = new Array<DP,4>(no_cylinder_shells+1,
				no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);				
		
		cylinder_ring_to_ring_Elsasser = new Array<DP,4>(no_cylinder_shells+1,
				no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);		
		
		cylinder_ring_to_ring_SF= new Array<DP,4>(no_cylinder_shells+1,
				no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);				
		
		*cylinder_ring_to_ring_self = 0.0;
		*cylinder_ring_to_ring_VF = 0.0;				
		*cylinder_ring_to_ring_Elsasser = 0.0;		
		*cylinder_ring_to_ring_SF = 0.0;
		
		forceV_cylinder_ring  = new Array<DP,2>(no_cylinder_shells+1,no_cylinder_kpll_slabs+1);					
		forceSF_cylinder_ring = new Array<DP,2>(no_cylinder_shells+1,no_cylinder_kpll_slabs+1);	
		
		*forceV_cylinder_ring = 0.0;					
		*forceSF_cylinder_ring = 0.0;
		
		energy_tr_cylinder_ring_B0 = new Array<DP,2>(no_cylinder_shells+1,
													no_cylinder_kpll_slabs+1);	
		*energy_tr_cylinder_ring_B0 = 0.0;
		
		
		if (ET_real_imag_switch == 1)
		{
			cylinder_ring_to_ring_self_real = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);
			
			cylinder_ring_to_ring_VF_real = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);	
			
			cylinder_ring_to_ring_Elsasser_real = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);			
			
			cylinder_ring_to_ring_SF_real = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);
			
			*cylinder_ring_to_ring_self_real = 0.0;
			*cylinder_ring_to_ring_VF_real = 0.0;	
			*cylinder_ring_to_ring_Elsasser_real = 0.0;			
			*cylinder_ring_to_ring_SF_real = 0.0;
			
			
			
			cylinder_ring_to_ring_self_imag = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);
			
			cylinder_ring_to_ring_VF_imag = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);	
			
			cylinder_ring_to_ring_Elsasser_imag = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);			
			
			cylinder_ring_to_ring_SF_imag = new Array<DP,4>(no_cylinder_shells+1,
					no_cylinder_kpll_slabs+1, no_cylinder_shells+1, no_cylinder_kpll_slabs+1);	
			
			*cylinder_ring_to_ring_self_imag = 0.0;
			*cylinder_ring_to_ring_VF_imag = 0.0;	
			*cylinder_ring_to_ring_Elsasser_imag = 0.0;			
			*cylinder_ring_to_ring_SF_imag = 0.0;	
			
			
			forceV_cylinder_ring_real	= new Array<DP,2>(no_cylinder_shells+1,
														no_cylinder_kpll_slabs+1);				
			forceSF_cylinder_ring_real	= new Array<DP,2>(no_cylinder_shells+1,
														no_cylinder_kpll_slabs+1);		
			forceV_cylinder_ring_imag	= new Array<DP,2>(no_cylinder_shells+1,
														no_cylinder_kpll_slabs+1);			
			forceSF_cylinder_ring_imag	= new Array<DP,2>(no_cylinder_shells+1,
														no_cylinder_kpll_slabs+1);	
			
			*forceV_cylinder_ring_real = 0.0;				
			*forceSF_cylinder_ring_real = 0.0;		
			*forceV_cylinder_ring_imag = 0.0;			
			*forceSF_cylinder_ring_imag = 0.0;	
		}
	}	// end of if (ET_anisotropic_cylinder_switch == 1)
	
}


//***************************   End  of EnergyTr.cc   *****************************************





